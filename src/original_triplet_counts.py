"""Validate unchanged original RNA triplets before enabling a processing notebook."""
import gzip
import hashlib
import io
from contextlib import ExitStack
import numpy as np
import pandas as pd
from scipy.io import mmread
from threadpoolctl import threadpool_limits


def explicit_rna_mask(genes, non_rna_feature_ids):
    """Select exact reviewed non-RNA IDs when a source omits modality labels."""
    ids = list(non_rna_feature_ids)
    if genes.shape[1] != 2 or not ids or len(set(ids)) != len(ids) or not set(ids).issubset(set(genes[0])):
        raise ValueError('Explicit non-RNA IDs require a two-column table and exact present unique IDs')
    return ~genes[0].isin(ids).to_numpy()


def audit(paths,seen,allowed_non_rna_types=(),allowed_feature_columns=(2,3),allow_plaintext_gz_members=False,exclude_ambiguous_barcodes=False,non_rna_feature_ids=()):
    encodings={}
    def table_bytes(role):
        raw=paths[role].read_bytes()
        compressed=raw.startswith(b'\x1f\x8b')
        if not compressed and not allow_plaintext_gz_members:
            raise ValueError('Plain-text member requires explicit source encoding review')
        encodings[role]='gzip' if compressed else 'reviewed_plaintext_despite_gz_suffix'
        return gzip.decompress(raw) if compressed else raw
    raw=table_bytes('genes_member')
    genes=pd.read_csv(io.BytesIO(raw),sep='\t',header=None,dtype=str,keep_default_na=False)
    bars=pd.read_csv(io.BytesIO(table_bytes('barcodes_member')),sep='\t',header=None,dtype=str,keep_default_na=False)[0]
    if not isinstance(exclude_ambiguous_barcodes,bool):
        raise ValueError('Explicit boolean ambiguous-barcode policy required')
    if not set(allowed_feature_columns).issubset({2,3,6}) or genes.shape[1] not in allowed_feature_columns or not genes[0].is_unique or genes[0].str.strip().eq('').any() or (not bars.is_unique and not exclude_ambiguous_barcodes) or bars.str.strip().eq('').any():
        raise ValueError('Invalid original gene/cell coordinates')
    rna=genes[2].eq('Gene Expression').to_numpy() if genes.shape[1]>=3 else np.ones(len(genes),dtype=bool)
    if non_rna_feature_ids:
        rna=explicit_rna_mask(genes,non_rna_feature_ids)
    non_rna_types=set(genes.loc[~rna,2]) if genes.shape[1]>=3 else set()
    if not rna.any() or not non_rna_types.issubset(set(allowed_non_rna_types)):
        raise ValueError('Mixed modalities require a separately reviewed RNA selection')
    with ExitStack() as stack,threadpool_limits(8):
        handle=stack.enter_context(paths['matrix_member'].open('rb'))
        compressed=handle.read(2)==b'\x1f\x8b';handle.seek(0)
        if not compressed and not allow_plaintext_gz_members:
            raise ValueError('Plain-text matrix requires explicit source encoding review')
        encodings['matrix_member']='gzip' if compressed else 'reviewed_plaintext_despite_gz_suffix'
        stream=stack.enter_context(gzip.GzipFile(fileobj=handle)) if compressed else handle
        x=mmread(stream,spmatrix=True).tocsc()
    x.check_format(full_check=True)
    if x.shape!=(len(genes),len(bars)) or not np.isfinite(x.data).all() or (x.data<0).any() or not np.equal(x.data,np.floor(x.data)).all():
        raise ValueError('Invalid original nonnegative integer RNA counts')
    x.sum_duplicates(); x.sort_indices(); x.eliminate_zeros()
    source_sha=hashlib.sha256(raw).hexdigest()
    source_features=len(genes)
    if non_rna_types or non_rna_feature_ids:
        lines=raw.splitlines(keepends=True)
        if len(lines)!=len(genes): raise ValueError('Ambiguous source feature-table rows')
        raw=b''.join(line for line,keep in zip(lines,rna) if keep)
        genes=genes.loc[rna].reset_index(drop=True)
        x=x[rna,:].tocsc(); x.sort_indices()
    sha=hashlib.sha256(raw).hexdigest()
    source_barcode_rows=len(bars)
    source_total_RNA_counts=int(x.sum())
    ambiguous=bars.duplicated(keep=False)
    excluded_ids=sorted(set(bars.loc[ambiguous])) if exclude_ambiguous_barcodes else []
    excluded_ambiguous_RNA_counts=int(x[:,ambiguous.to_numpy()].sum()) if exclude_ambiguous_barcodes else 0
    if exclude_ambiguous_barcodes:
        # Validate every original value/coordinate above before excluding all
        # columns carrying a repeated identifier. Never keep an arbitrary copy.
        keep=~ambiguous.to_numpy()
        x=x[:,keep].tocsc();bars=bars.loc[keep].reset_index(drop=True)
        if not len(bars):raise ValueError('No unambiguous original barcode columns')
    # Avoid iterating millions of empty droplet whitelist rows in Python.
    for i in np.flatnonzero(np.diff(x.indptr)>=20):
        barcode=bars.iloc[i]
        start,end=x.indptr[i:i+2]
        if end-start<20: continue
        digest=hashlib.sha256(x.indices[start:end].astype('<i8').tobytes()+x.data[start:end].astype('<i8').tobytes()).hexdigest()
        key=(barcode,sha,digest)
        if key in seen: raise ValueError('Exact original cell/count duplication across captures')
        seen.add(key)
    return dict(cells=len(bars),genes=len(genes),nonzeros=int(x.nnz),gene_coordinate_sha256=sha,
        integer_nonnegative_counts=True,original_gene_symbols_preserved=True,
        source_feature_table_sha256=source_sha,source_features=source_features,
        source_barcode_rows=source_barcode_rows,RNA_positive_barcode_rows=int(np.count_nonzero(np.diff(x.indptr))),
        excluded_ambiguous_barcode_rows=source_barcode_rows-len(bars),excluded_ambiguous_barcode_ids=excluded_ids,
        source_total_RNA_counts=source_total_RNA_counts,excluded_ambiguous_RNA_counts=excluded_ambiguous_RNA_counts,
        barcode_scope='Source matrix barcode rows; not automatically called or viable cells',
        excluded_non_rna_features=source_features-len(genes),reviewed_non_rna_types=sorted(non_rna_types),
        explicit_non_rna_feature_ids=list(non_rna_feature_ids),source_member_encodings=encodings)
