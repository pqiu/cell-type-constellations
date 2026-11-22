"""Pin original pure-RNA H5 axes and QC before count-preserving conversion."""
import hashlib

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

from scripts.download_pbmc_asset import ROOT
from src.lyme_vireo import audit_h5


def digest(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def audit_RNA_h5(path, accession, cohort, seen, expected_genome, allowed_non_RNA_types=(), source_gene_id_mode='ensembl'):
    """Validate all original values; retain source QC for the fixed RNA gate."""
    folder = ROOT / 'references/datasets' / (cohort + '_original_RNA_QC')
    folder.mkdir(exist_ok=True)
    with h5py.File(path, 'r') as h:
        matrix = h['matrix']; features = matrix['features']
        types = features['feature_type'].asstr()[()]
        rna = types == 'Gene Expression'
        if not rna.any() or not (set(types) - {'Gene Expression'}).issubset(set(allowed_non_RNA_types)):
            raise ValueError('Non-RNA source types require an explicit source-specific review')
        ids = features['id'].asstr()[()][rna]; symbols = features['name'].asstr()[()][rna]
        if not pd.Index(ids).is_unique:
            raise ValueError('Original unique RNA coordinates required')
        if source_gene_id_mode=='ensembl':
            if not pd.Series(ids).str.fullmatch(r'ENSG[0-9]+(?:\.[0-9]+)?').all():
                raise ValueError('Original unique Ensembl RNA coordinates required')
        elif source_gene_id_mode=='native_symbols':
            if not np.array_equal(ids,symbols) or any(not str(v).strip() for v in ids):
                raise ValueError('Explicit native symbol IDs must equal original names and be nonempty')
        else:
            raise ValueError('Unknown source RNA identifier mode')
        if set(features['genome'].asstr()[()][rna]) != {expected_genome}:
            raise ValueError('Native genome label differs from the source-specific review')
        mito = np.array([s.startswith('MT-') for s in symbols])
        if int(mito.sum()) != 13:
            raise ValueError('Primary mitochondrial RNA coordinates need separate review')
        model = set((ROOT / 'models/model_v1.1/gene_order.tsv').read_text().splitlines())
        overlap = len(model & set(symbols))
        if overlap < 5000:
            raise ValueError('Insufficient native SCimilarity symbol compatibility')
        bars = matrix['barcodes'].asstr()[()]
        if not pd.Index(bars).is_unique:
            raise ValueError('Ambiguous original RNA barcodes')
        values = matrix['data'][()]
        if not np.isfinite(values).all() or (values < 0).any() or not np.equal(values, np.floor(values)).all():
            raise ValueError('Invalid original RNA count values')
        counts = sparse.csc_matrix((values,matrix['indices'][()],matrix['indptr'][()]),
            shape=tuple(matrix['shape'][()])).T.tocsr()
    counts.check_format(full_check=True)
    counts.sum_duplicates(); counts.sort_indices(); counts.eliminate_zeros()
    if counts.shape != (len(bars),len(types)):
        raise ValueError('Native count and coordinate dimensions disagree')
    counts = counts[:,rna].tocsr()
    totals = np.asarray(counts.sum(axis=1,dtype=np.int64)).ravel()
    keep = totals >= 500
    selected = counts[keep]
    table = pd.DataFrame(dict(source_barcode=bars[keep],RNA_total=totals[keep],
        RNA_detected_genes=np.diff(selected.indptr),
        mito_RNA_total=np.asarray(selected[:,mito].sum(axis=1,dtype=np.int64)).ravel()))
    review = audit_h5(path,pd.DataFrame(index=bars[keep]),seen)
    if review['shape'] != list(counts.shape):
        raise ValueError('Independent native RNA audit disagrees')
    qc_path = folder / (accession + '.csv.gz')
    if qc_path.exists(): pd.testing.assert_frame_equal(pd.read_csv(qc_path),table,check_dtype=False)
    else: table.to_csv(qc_path,index=False,compression='gzip')
    axis = pd.DataFrame(dict(source_gene_id=ids,source_gene_symbols=symbols,annotation_gene_symbols=symbols))
    axis_path = folder / (review['gene_coordinate_sha256'] + '_genes.csv.gz')
    if axis_path.exists(): pd.testing.assert_frame_equal(pd.read_csv(axis_path,keep_default_na=False),axis)
    else: axis.to_csv(axis_path,index=False,compression='gzip')
    review.update(source_barcode_rows=len(bars),cells=len(bars),genes=len(ids),source_gene_id_mode=source_gene_id_mode,
        source_feature_types=pd.Series(types).value_counts().to_dict(),source_genome_label=expected_genome,
        excluded_non_RNA_features=int((~rna).sum()),reviewed_non_RNA_types=sorted(set(types[~rna])),
        annotation_symbol_overlap=overlap,
        retained_native_coordinates=dict(path=axis_path.relative_to(ROOT).as_posix(),sha256=digest(axis_path)),
        preconversion_RNA_QC_review=dict(source_barcode_rows=len(bars),minimum_RNA_counts=500,
            preconversion_selected_barcodes=int(keep.sum()),full_source_RNA_counts=int(totals.sum()),
            selected_RNA_counts=int(totals[keep].sum()),all_original_values_and_coordinates_checked=True,
            metadata=dict(path=qc_path.relative_to(ROOT).as_posix(),sha256=digest(qc_path))))
    return review


def audit_pure_RNA(path, accession, cohort, seen, expected_genome):
    """Strict convenience entry point; additional modalities are never implicit."""
    return audit_RNA_h5(path,accession,cohort,seen,expected_genome)
