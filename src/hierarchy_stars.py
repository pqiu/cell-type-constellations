"""Construct ontology-parent centroids from exclusive annotation labels.

An exact parent mean is the cell-count-weighted mean of exclusive direct-label
means, including direct-label groups too small to have had an exported star.
Each source cell is counted once per parent even when several ontology paths
connect its label to that parent.
"""
from __future__ import annotations

import numpy as np
from scipy import sparse


def log_normalized_sum(counts, gene_columns):
    """Sum per-cell log1p(TP10K) on selected genes, with full-source library sizes."""
    matrix = sparse.csr_matrix(counts, dtype=np.float64)
    libraries = np.asarray(matrix.sum(axis=1)).ravel()
    if not np.all(np.isfinite(libraries)) or np.any(libraries <= 0):
        raise ValueError("All included cells need positive finite full-source libraries")
    selected = matrix[:, np.asarray(gene_columns, dtype=np.int64)].tocsr()
    normalized = selected.multiply((10000 / libraries)[:, None]).tocsr()
    normalized.data = np.log1p(normalized.data)
    return np.asarray(normalized.sum(axis=0)).ravel()


def pooled_parent_means(direct_sums, direct_counts, incidence, min_cells=20):
    """Return parent means and support from exclusive direct-label summaries.

    ``incidence`` has one row per direct SCimilarity label and one column per
    requested parent. Multiple graph paths must already be deduplicated.
    """
    sums = np.asarray(direct_sums, dtype=np.float64)
    counts = np.asarray(direct_counts, dtype=np.int64)
    links = sparse.csr_matrix(incidence, dtype=np.int64)
    if sums.ndim != 2 or len(counts) != len(sums) or links.shape[0] != len(counts):
        raise ValueError("Incompatible direct-label inputs")
    if np.any(counts < 0) or np.any((links.data != 1)):
        raise ValueError("Negative counts or nondeduplicated membership")
    parent_counts = np.asarray(counts @ links).ravel()
    parent_sums = np.asarray(links.T @ sums)
    parent_means = np.full(parent_sums.shape, np.nan, dtype=np.float64)
    supported = parent_counts >= min_cells
    parent_means[supported] = parent_sums[supported] / parent_counts[supported, None]
    return parent_means, parent_counts
