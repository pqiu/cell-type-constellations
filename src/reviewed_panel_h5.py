"""Read reviewed source RNA panels without treating excluded genes as zeros."""
import hashlib
import json
import re
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

from scripts.download_pbmc_asset import ROOT, LIMIT, footprint, download


def digest(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def read_counts(path, minimum_RNA_counts=500):
    if minimum_RNA_counts != 500:
        raise ValueError('Explicitly reviewed500RNAcount gate required')
    with h5py.File(path) as h:
        matrix = h['matrix']; features = matrix['features']
        ids = pd.Index(features['id'].asstr()[:]); symbols = features['name'].asstr()[:]
        types = features['feature_type'].asstr()[:]; rna = types == 'Gene Expression'
        if not ids.is_unique or not rna.any() or not set(types) <= {'Gene Expression', 'Antibody Capture'}:
            raise ValueError('Unexpected original feature coordinates/types')
        if set(features['genome'].asstr()[:][rna]) != {'GRCh38'}:
            raise ValueError('Reviewed human GRCh38 source required')
        canonical = np.asarray([bool(re.fullmatch(r'ENSG[0-9]+(?:\.[0-9]+)?', gene)) for gene in ids]) & rna
        if any(not str(gene).startswith('DEPRECATED_') for gene in ids[rna & ~canonical]):
            raise ValueError('Unknown noncanonical RNA feature identity')
        sets = list(features.get('target_sets', {}))
        if sets:
            if len(sets) != 1:
                raise ValueError('Exactly one reviewed source target set required')
            panel, = sets
            values = features['target_sets'][panel].asstr()[:]
            if not all(re.fullmatch(r'[0-9]+', value) for value in values):
                raise ValueError('Source target feature indices required')
            indices = values.astype(np.int64)
            if not len(indices) or len(set(indices)) != len(indices) or indices.min() < 0 or indices.max() >= len(ids):
                raise ValueError('Invalid source target indices')
            keep_genes = np.zeros(len(ids), dtype=bool); keep_genes[indices] = True
            if (keep_genes & ~canonical).any():
                raise ValueError('Source-approved targets include noncanonical/non-RNA coordinates')
            mode = 'source_approved_probe_target_set'
        else:
            panel = 'no_source_target_set'; keep_genes = canonical
            if not np.array_equal(keep_genes, rna):
                raise ValueError('Unrestricted source has unreviewed noncanonical RNA rows')
            mode = 'all_canonical_source_RNA_genes'
        barcodes = pd.Index(matrix['barcodes'].asstr()[:])
        counts = sparse.csc_matrix((matrix['data'][:], matrix['indices'][:], matrix['indptr'][:]), shape=tuple(matrix['shape'][:]))
        counts.check_format(full_check=True)
        if not counts.has_canonical_format or counts.shape != (len(ids), len(barcodes)) or not barcodes.is_unique:
            raise ValueError('Invalid original sparse/barcode coordinates')
        if not np.isfinite(counts.data).all() or (counts.data < 0).any() or not np.equal(counts.data, np.floor(counts.data)).all():
            raise ValueError('All original values, including excluded features, must be nonnegative integer counts')
        counts.eliminate_zeros()
        # Keep original selected counts and original feature order; never alias
        # deprecated probes back to canonical genes or fabricate excluded zeros.
        panel_counts = counts[keep_genes].tocsc()
        totals = np.asarray(panel_counts.sum(axis=0, dtype=np.int64)).ravel()
        support = np.diff(panel_counts.indptr)
        mito_genes = np.char.startswith(symbols[keep_genes].astype(str), 'MT-')
        if not mito_genes.any():
            raise ValueError('Mitochondrial gene coverage needs separate QC review')
        mito = np.asarray(panel_counts[mito_genes].sum(axis=0, dtype=np.int64)).ravel()
        fraction = np.divide(mito, totals, out=np.zeros_like(totals, dtype=float), where=totals > 0)
        chosen = totals >= minimum_RNA_counts
        full_RNA = np.asarray(counts[rna].sum(axis=0, dtype=np.int64)).ravel()
        full_features = np.asarray(counts.sum(axis=0, dtype=np.int64)).ravel()
        obs = pd.DataFrame(dict(source_original_barcode_index=np.flatnonzero(chosen),
            source_geometry_RNA_counts=totals[chosen], source_geometry_detected_genes=support[chosen],
            source_geometry_mito_RNA_counts=mito[chosen], source_geometry_mito_fraction=fraction[chosen],
            source_original_RNA_counts=full_RNA[chosen], source_original_all_feature_counts=full_features[chosen]), index=barcodes[chosen])
        var = pd.DataFrame(dict(gene_symbols=symbols[keep_genes], feature_type='Gene Expression',
            source_original_feature_index=np.flatnonzero(keep_genes)), index=ids[keep_genes])
        metadata = dict(source_gene_policy=mode, source_target_set=panel,
            original_features=len(ids), original_RNA_features=int(rna.sum()), geometry_genes=int(keep_genes.sum()),
            excluded_non_RNA_features=int((~rna).sum()), excluded_RNA_features=int((rna & ~keep_genes).sum()),
            original_barcode_rows=len(barcodes), stage_RNA_barcodes=int(chosen.sum()),
            original_all_feature_UMI=int(full_features.sum()), original_RNA_UMI=int(full_RNA.sum()),
            geometry_RNA_UMI=int(totals.sum()), excluded_RNA_UMI=int((full_RNA - totals).sum()),
            geometry_mito_genes=int(mito_genes.sum()),
            gene_coordinate_sha256=hashlib.sha256('\n'.join(ids[keep_genes]).encode()).hexdigest(),
            provisional_project_QC_cells=int(((totals >= 500) & (support >= 200) & (fraction <= .2)).sum()))
        for key, field in [('source_library_ids', 'library_ids'), ('source_chemistry', 'chemistry_description'), ('source_software', 'software_version')]:
            value = h.attrs[field]
            if isinstance(value, (bytes, str)):
                metadata[key] = value.decode() if isinstance(value, bytes) else value
            else:
                metadata[key] = [v.decode() if isinstance(v, bytes) else str(v) for v in value]
        return panel_counts[:, chosen].T.tocsr(), var, obs, metadata


def prepare(cfg, destination):
    asset = cfg['source_asset']; name = asset['filename']
    if Path(name).name != name:
        raise ValueError('Original source basename required')
    path = download(asset['url'], ROOT / 'data/downloads' / (cfg['study_id'] + '_assets') / name, asset['expected_bytes'])
    if digest(path) != asset['sha256']:
        raise ValueError('Original count source differs from pinned native audit')
    counts, var, obs, review = read_counts(path)
    source = cfg['source_metadata']['original_count_review']
    checks = dict(original_RNA_features='original_RNA_features', geometry_genes='target_set_features',
        original_barcode_rows='original_barcode_rows', stage_RNA_barcodes='target_set_500RNA_barcodes',
        original_all_feature_UMI='original_all_feature_UMI', original_RNA_UMI='original_RNA_UMI',
        geometry_RNA_UMI='target_set_UMI', excluded_RNA_UMI='excluded_if_target_filtered_UMI',
        provisional_project_QC_cells='target_set_provisional_project_QC_cells', gene_coordinate_sha256='target_set_gene_coordinate_sha256')
    for target, original in checks.items():
        if review[target] != source[original]:
            raise ValueError('Native audit mismatch: ' + target)
    if review['source_target_set'] != source['target_set_name'] or review['source_library_ids'] != [source['source_internal_library']]:
        raise ValueError('Source library or approved gene set differs from native review')
    if footprint() + 2 * (counts.data.nbytes + counts.indices.nbytes + counts.indptr.nbytes) + 2_000_000_000 >= LIMIT:
        raise RuntimeError('RNA panel conversion exceeds project disk budget')
    ad.AnnData(X=counts, obs=obs, var=var).write_h5ad(destination, compression='gzip')
    manifest = json.loads(path.with_suffix(path.suffix + '.manifest.json').read_text())
    return dict(kind='reviewed_source_RNA_gene_panel', asset=manifest, conversion=review,
        geometry_denominator='Retained approved RNA gene set per cell; antibody, excluded and deprecated features excluded. Absent genes remain unavailable outside saved gene axis.',
        QC_scope='Project500RNA/200genes/20percentMT applied to retained RNA set; MT coverage recorded explicitly; not a claim of reproducing author cell calls.'), [path]
