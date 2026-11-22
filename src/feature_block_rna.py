"""Select explicitly typed RNA rows from a canonical original feature block."""
import numpy as np
from scipy import sparse


def select_rna_rows(matrix, labels, feature_types):
    matrix=sparse.csr_matrix(matrix)
    if matrix.shape[0]!=len(labels) or len(labels)!=len(feature_types):
        raise ValueError('Feature labels and types must cover every original block row')
    if not matrix.has_canonical_format:
        raise ValueError('Original sparse feature rows must be canonical')
    if not np.isfinite(matrix.data).all() or (matrix.data<0).any() or not np.equal(matrix.data,np.floor(matrix.data)).all():
        raise ValueError('Original feature values must be nonnegative integer counts')
    if set(feature_types)-{'Gene Expression','Antibody Capture'}:
        raise ValueError('Unreviewed source feature modality')
    selected=np.asarray(feature_types)=='Gene Expression'
    return matrix[selected].tocsr(),np.asarray(labels,dtype=str)[selected]
