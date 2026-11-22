"""Bound allocated star chunks rather than the unallocated NaN-filled shape."""


def output_gene_indices(symbols, panel):
    """Return source-axis indices for an ordered, unique output-gene panel.

    This only controls persisted star vectors.  Annotation, QC, and library-size
    normalization continue to use the complete source count axis.
    """
    from collections import Counter
    import numpy as np

    source = [str(value) for value in symbols]
    requested = [str(value) for value in panel]
    if not requested or len(requested) != len(set(requested)):
        raise ValueError('Output gene panel must be nonempty and unique')
    frequency = Counter(source)
    ambiguous = [gene for gene in requested if frequency[gene] != 1]
    if ambiguous:
        raise ValueError(
            'Every output-panel gene must occur exactly once on the source axis: '
            + ','.join(ambiguous[:10])
        )
    lookup = {gene: index for index, gene in enumerate(source)}
    return np.asarray([lookup[gene] for gene in requested], dtype=np.int64)


def export_budget(labels,views,min_cells,n_genes,n_classes,obs_bytes=0,availability_bytes=0,view_count=None,array_count=10):
    import pandas as pd
    counts=pd.DataFrame(dict(view=list(views),label=list(labels))).groupby(['view','label'],observed=True).size()
    supported=int(counts.ge(min_cells).sum())
    observed=len(counts)
    n_views=len(set(views)) if view_count is None else int(view_count)
    if n_views<len(set(views)):raise ValueError('Budget view count omits retained source views')
    if not isinstance(array_count, int) or array_count < 1 or array_count > 10:
        raise ValueError('array_count must be an integer from 1 through 10')
    # Each float64 array allocates one full gene-row chunk only when supported.
    # Allow4KiB perchunk forcompression/index overhead, plus metadata/CSV space.
    gene_support_rows=observed if availability_bytes else 0
    array_bytes=(supported*array_count+gene_support_rows)*(n_genes*8+4096)
    metadata_bytes=2*int(obs_bytes)+n_views*n_classes*2048+100_000_000
    return dict(supported_stars=supported,observed_types=observed,array_count=array_count,
                expected_bytes=int(array_bytes+metadata_bytes+availability_bytes))
