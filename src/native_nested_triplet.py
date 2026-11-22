"""Audit unchanged count triplets contained in public sample archives."""
import gzip
import hashlib
import io
import tarfile
import numpy as np
import pandas as pd
from scipy.io import mmread
from threadpoolctl import threadpool_limits


def audit(path, seen, allowed_non_rna_types=()):
    members = {}; data = {}
    with tarfile.open(path) as archive:
        for role,endings in [('genes_member',('features.tsv.gz','genes.tsv.gz','features.tsv','genes.tsv')),
                             ('barcodes_member',('barcodes.tsv.gz','barcodes.tsv')),
                             ('matrix_member',('matrix.mtx.gz','matrix.mtx'))]:
            hits = [m for m in archive.getmembers() if m.isfile() and m.name.endswith(endings)]
            if len(hits) != 1 or hits[0].size > 5_000_000_000:
                raise ValueError('Ambiguous or oversized native nested count member')
            member = hits[0];members[role] = member.name
            raw = archive.extractfile(member).read()
            data[role] = gzip.decompress(raw) if raw.startswith(b'\x1f\x8b') else raw
    genes = pd.read_csv(io.BytesIO(data['genes_member']),sep='\t',header=None,dtype=str,keep_default_na=False)
    bars = pd.read_csv(io.BytesIO(data['barcodes_member']),sep='\t',header=None,dtype=str,keep_default_na=False)
    if (genes.shape[1] not in [2,3] or bars.shape[1] != 1 or not genes[0].is_unique
            or not bars[0].is_unique or genes[0].str.strip().eq('').any() or bars[0].str.strip().eq('').any()):
        raise ValueError('Invalid native gene/cell coordinates')
    rna = genes[2].eq('Gene Expression').to_numpy() if genes.shape[1] == 3 else np.ones(len(genes),dtype=bool)
    types = set(genes.loc[~rna,2]) if genes.shape[1] == 3 else set()
    if not rna.any() or not types.issubset(set(allowed_non_rna_types)):
        raise ValueError('Non-RNA features require explicit source review')
    with threadpool_limits(8):
        x = mmread(io.BytesIO(data['matrix_member']),spmatrix=True).tocsc()
    x.check_format(full_check=True)
    if (x.shape != (len(genes),len(bars)) or not np.isfinite(x.data).all()
            or (x.data < 0).any() or not np.equal(x.data,np.floor(x.data)).all()):
        raise ValueError('Invalid original nonnegative integer counts')
    x.sum_duplicates();x.sort_indices();x.eliminate_zeros()
    native_features = len(genes)
    raw = data['genes_member']
    if types:
        lines = raw.splitlines(keepends=True)
        if len(lines) != len(genes):
            raise ValueError('Ambiguous original feature-table rows')
        raw = b''.join(line for line,keep in zip(lines,rna) if keep)
        genes = genes.loc[rna].reset_index(drop=True)
        x = x[rna,:].tocsc();x.sort_indices()
    sha = hashlib.sha256(raw).hexdigest()
    for i in np.flatnonzero(np.diff(x.indptr) >= 20):
        start,end = x.indptr[i:i+2]
        digest = hashlib.sha256(x.indices[start:end].astype('<i8').tobytes()+x.data[start:end].astype('<i8').tobytes()).hexdigest()
        key = (bars.iloc[i,0],sha,digest)
        if key in seen:
            raise ValueError('Exact native cell/count duplication across source archives')
        seen.add(key)
    review = dict(cells=len(bars),genes=len(genes),nonzeros=int(x.nnz),integer_nonnegative_counts=True,
        source_barcode_rows=len(bars),gene_coordinate_sha256=sha,
        excluded_non_RNA_features=native_features-len(genes),reviewed_non_RNA_types=sorted(types),
        member_sha256={role:hashlib.sha256(raw).hexdigest() for role,raw in data.items()},
        source_mitochondrial_features=int(genes[1].str.startswith('MT-').sum()))
    return members,review,genes
