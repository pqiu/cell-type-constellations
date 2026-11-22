"""Select explicitly reviewed human gene coordinates from a mixed count assay."""
import hashlib
import json

import numpy as np
import pandas as pd
from scipy import sparse


def coordinate_sha(ids, symbols):
    return hashlib.sha256(json.dumps([list(ids), list(symbols)], separators=(',', ':')).encode()).hexdigest()


def feature_sha(genes):
    return hashlib.sha256(json.dumps(genes.astype(str).values.tolist(), separators=(',', ':')).encode()).hexdigest()


def select(counts, genes, review):
    """Validate every source value before retaining pinned native human genes."""
    counts = sparse.csr_matrix(counts)
    counts.check_format(full_check=True)
    if genes.shape[1] != 3 or counts.shape[1] != len(genes):
        raise ValueError('Three-column native source feature axis required')
    values = counts.data
    if not np.isfinite(values).all() or (values < 0).any() or not np.equal(values, np.floor(values)).all():
        raise ValueError('Invalid original counts, including excluded features')
    counts.sum_duplicates(); counts.sort_indices(); counts.eliminate_zeros()
    if feature_sha(genes) != review['source_feature_axis_sha256']:
        raise ValueError('Pinned original feature axis changed')
    types = set(genes[2])
    if types != set(review['source_feature_types']) or not types.issubset({'Gene Expression', 'Antibody Capture'}):
        raise ValueError('Unreviewed original feature modalities')
    rna = genes[2].eq('Gene Expression').to_numpy()
    human = rna & genes[0].str.fullmatch(r'ENSG[0-9]+(?:\.[0-9]+)?').to_numpy()
    if not human.any() or not genes.loc[human, 0].is_unique:
        raise ValueError('Original unique human Ensembl coordinates required')
    excluded = genes.loc[rna & ~human, 0].tolist()
    if excluded != review['excluded_RNA_feature_ids']:
        raise ValueError('Unreviewed nonhuman RNA feature IDs')
    native = genes.loc[human].reset_index(drop=True)
    if len(native) != review['human_gene_rows'] or coordinate_sha(native[0], native[1]) != review['human_gene_coordinate_sha256']:
        raise ValueError('Pinned human gene coordinates changed')
    libraries = pd.DataFrame(dict(
        source_all_RNA_counts=np.asarray(counts[:, rna].sum(axis=1, dtype=np.int64)).ravel(),
        source_excluded_nonhuman_RNA_counts=np.asarray(counts[:, rna & ~human].sum(axis=1, dtype=np.int64)).ravel(),
        source_antibody_counts=np.asarray(counts[:, ~rna].sum(axis=1, dtype=np.int64)).ravel()))
    return counts[:, human].tocsr(), native, libraries
