"""Read materialized selected-gene stars and native-gene star recipes.

Native parent profiles are reconstructed from the original supported direct
label means and saved sums for direct labels with 1–19 cells. This preserves
every computable ancestor profile without storing redundant dense native
vectors for overlapping parents. Reconstruction agrees with materialized
selected-gene profiles to their float32 precision.
"""
from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import pandas as pd


def reconstruct_native_parent(root: Path, dataset_id: str, sample_id: str,
                              parent_id: str) -> tuple[int, np.ndarray]:
    root=Path(root)
    curated=root/'results/star_maps'/dataset_id
    rebuilt=root/'results/hierarchy_all_star_maps'/dataset_id
    stars=pd.read_csv(rebuilt/'stars.csv',usecols=['sample_id','parent_id','n_cells','status'])
    target=stars.loc[stars.sample_id.eq(sample_id)&stars.parent_id.eq(parent_id)]
    if len(target)!=1:
        raise KeyError((dataset_id,sample_id,parent_id))
    n=int(target.iloc[0].n_cells)
    if n==0:
        raise ValueError('Parent star is absent from this map')
    membership=pd.read_csv(root/'references/vocabulary/aggregation_membership.tsv',sep='\t')
    labels=set(membership.loc[membership.aggregate_id.eq(parent_id),'original_label'])
    if not labels:
        raise KeyError(parent_id)
    direct=pd.read_csv(curated/'stars.csv',usecols=['sample_id','label','n_cells','status'])
    direct['source_row']=np.arange(len(direct),dtype=np.int64)
    direct=direct.loc[direct.sample_id.eq(sample_id)&direct.label.isin(labels)&direct.n_cells.gt(0)]
    if direct.n_cells.sum()!=n:
        raise AssertionError('Direct support does not match parent support')
    basis=pd.read_csv(rebuilt/'under20_direct_rows.csv')
    basis=basis.loc[basis.sample_id.eq(sample_id)].set_index('label')
    direct_means = rebuilt/'native_direct_mean_basis.h5'
    if not direct_means.exists():
        direct_means = curated/'expression.h5'
    with h5py.File(direct_means,'r') as old, h5py.File(
            rebuilt/'under20_direct_sums_native.h5','r') as small:
        source=old['norm10k_log1p_mean']
        sums=small['normalized_log_sum']
        total=np.zeros(source.shape[1],dtype=np.float64)
        for row in direct.itertuples(index=False):
            if row.status=='supported':
                total+=np.asarray(source[row.source_row],dtype=np.float64)*row.n_cells
            else:
                entry=basis.loc[row.label]
                if int(entry.n_cells)!=row.n_cells:
                    raise AssertionError('Small direct-label support changed')
                total+=np.asarray(sums[int(entry.basis_row)],dtype=np.float64)
    return n,(total/n).astype(np.float32)
