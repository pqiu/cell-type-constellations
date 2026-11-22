"""Keep a centroid undefined if any contributing cell lacks a source gene."""
import numpy as np
from scipy import sparse

def validate_availability(groups, availability, n_genes):
    groups = np.asarray(groups)
    availability = np.asarray(availability)
    if groups.ndim != 1 or groups.dtype.kind not in 'iu':
        raise ValueError('Integer source-library indices required')
    if availability.ndim != 2 or availability.shape[1] != n_genes or availability.dtype != np.dtype(bool):
        raise ValueError('Boolean source-library by gene availability required')
    if groups.size and (groups.min() < 0 or groups.max() >= availability.shape[0]):
        raise ValueError('Source-library index outside availability table')
    return groups, availability

def observed_gene_support(groups, availability, n_genes):
    groups, availability = validate_availability(groups, availability, n_genes)
    frequencies = np.bincount(groups, minlength=len(availability))
    return frequencies @ availability.astype(np.int64)

def require_observed_counts(matrix, groups, availability):
    x = sparse.csr_matrix(matrix)
    groups, availability = validate_availability(groups, availability, x.shape[1])
    if len(groups) != x.shape[0]:
        raise ValueError('Cell/source-library dimensions differ')
    for library in np.unique(groups):
        if x[groups == library][:,~availability[library]].count_nonzero():
            raise ValueError('Nonzero counts in source-unavailable genes')

def mask_incomplete_summary(summary, groups, availability):
    """Retain exact all-cell statistics; never average only observed cells."""
    n_genes = len(summary[0])
    support = observed_gene_support(groups, availability, n_genes)
    n = len(groups)
    incomplete = support != n
    result = tuple(np.asarray(values,dtype=np.float64).copy() for values in summary)
    for values in result:
        if values.shape != (n_genes,):
            raise ValueError('Summary/gene dimensions differ')
        values[incomplete] = np.nan
    return result, support
