"""Shared source-partial calculation for within-map gene variance."""
import json
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from scripts.select_hierarchy_all_genes import ROOT, SOURCE, reconstruct_sample
from scripts.star_matrix_core import digest, symbol_groups

PARTIALS = ROOT / "results/hierarchy_global_reanalysis/final_132/pbmc_66_within_map_gene_selection_v1/source_partials"

def source_partial(
    dataset_id: str,
    rows: pd.DataFrame,
    genes: list[str],
    source_labels: list[str],
    parents: list[str],
    incidence,
) -> Path:
    PARTIALS.mkdir(parents=True, exist_ok=True)
    path = PARTIALS / f"{dataset_id}.npz"
    groups = symbol_groups(ROOT / "results/star_maps" / dataset_id / "genes.csv")
    present_genes = [gene for gene in genes if gene in groups]
    gene_position = {gene: i for i, gene in enumerate(genes)}
    gene_indices = np.asarray([gene_position[gene] for gene in present_genes], dtype=np.int32)
    expanded_rows = np.concatenate([groups[gene] for gene in present_genes])
    expanded_names = [f"__native_row_{i}" for i in range(len(expanded_rows))]
    expanded_axis = dict(zip(expanded_names, map(int, expanded_rows)))
    duplicate_counts = np.asarray([len(groups[gene]) for gene in present_genes], dtype=np.int64)
    duplicate_starts = np.r_[0, np.cumsum(duplicate_counts)[:-1]].astype(np.int64)

    expected_maps = int(rows["map_index"].nunique())
    if path.exists():
        with np.load(path, allow_pickle=False) as saved:
            if (
                int(saved["retained_stars"]) == len(rows)
                and int(saved["retained_maps"]) == expected_maps
                and np.array_equal(saved["gene_indices"], gene_indices)
            ):
                return path
        raise ValueError(f"Stale within-map partial: {path}")

    curated = ROOT / "results/star_maps" / dataset_id
    hierarchy = SOURCE / dataset_id
    direct = pd.read_csv(curated / "stars.csv", usecols=["sample_id", "label", "n_cells", "status"])
    direct["direct_row"] = np.arange(len(direct), dtype=np.int64)
    small = pd.read_csv(hierarchy / "under20_direct_rows.csv")
    old_path = curated / "expression.h5"
    if not old_path.exists():
        old_path = hierarchy / "native_direct_mean_basis.h5"
        receipt = json.loads((hierarchy / "retired_direct_label_expression.json").read_text())
        if digest(old_path) != receipt["basis_sha256"]:
            raise ValueError(f"Retired direct basis hash changed: {dataset_id}")

    parent_position = {name: i for i, name in enumerate(parents)}
    variance_sum = np.zeros(len(present_genes), dtype=np.float64)
    variance_sumsq = np.zeros(len(present_genes), dtype=np.float64)
    valid_maps = np.zeros(len(present_genes), dtype=np.int32)
    minimum_finite_stars = np.full(len(present_genes), np.iinfo(np.int16).max, dtype=np.int16)
    maximum_finite_stars = np.zeros(len(present_genes), dtype=np.int16)
    processed_maps = 0
    with (
        h5py.File(old_path) as old_h5,
        h5py.File(hierarchy / "under20_direct_sums_native.h5") as small_h5,
    ):
        for sample_id, sample_rows in rows.groupby("sample_id", sort=False):
            if sample_rows["map_index"].nunique() != 1:
                raise ValueError(f"A source sample maps to multiple maps: {dataset_id} {sample_id}")
            means, parent_counts = reconstruct_sample(
                dataset_id,
                sample_id,
                expanded_names,
                expanded_axis,
                source_labels,
                incidence,
                direct,
                small,
                old_h5,
                small_h5,
            )
            positions = np.asarray([parent_position[x] for x in sample_rows.parent_id], dtype=np.int64)
            if not np.array_equal(parent_counts[positions], sample_rows.n_cells.to_numpy(dtype=np.int64)):
                raise ValueError(f"Parent support changed: {dataset_id} {sample_id}")
            expanded_values = means[positions]
            finite_expanded = np.isfinite(expanded_values)
            collapsed = np.add.reduceat(
                np.where(finite_expanded, expanded_values, 0.0), duplicate_starts, axis=1
            )
            observed_duplicate_rows = np.add.reduceat(
                finite_expanded.astype(np.int16), duplicate_starts, axis=1
            )
            values = np.divide(
                collapsed,
                observed_duplicate_rows,
                out=np.full(collapsed.shape, np.nan, dtype=np.float64),
                where=observed_duplicate_rows > 0,
            )
            finite = np.isfinite(values)
            counts = finite.sum(axis=0, dtype=np.int64)
            usable = counts >= 2
            if usable.any():
                totals = np.where(finite, values, 0.0).sum(axis=0, dtype=np.float64)
                totals_sq = np.where(finite, np.square(values), 0.0).sum(axis=0, dtype=np.float64)
                centered_ss = totals_sq[usable] - np.square(totals[usable]) / counts[usable]
                tolerance = np.maximum(1.0, totals_sq[usable]) * 1e-12
                if np.any(centered_ss < -tolerance):
                    raise ValueError(f"Negative centered sum of squares: {dataset_id} {sample_id}")
                variances = np.maximum(centered_ss, 0.0) / counts[usable]
                variance_sum[usable] += variances
                variance_sumsq[usable] += np.square(variances)
                valid_maps[usable] += 1
                minimum_finite_stars[usable] = np.minimum(
                    minimum_finite_stars[usable], counts[usable].astype(np.int16)
                )
                maximum_finite_stars[usable] = np.maximum(
                    maximum_finite_stars[usable], counts[usable].astype(np.int16)
                )
            processed_maps += 1

    if processed_maps != expected_maps:
        raise ValueError(f"Map count changed for {dataset_id}: {processed_maps} != {expected_maps}")
    minimum_finite_stars[valid_maps == 0] = 0
    temporary = path.with_suffix(".building.npz")
    np.savez_compressed(
        temporary,
        gene_indices=gene_indices,
        variance_sum=variance_sum,
        variance_sumsq=variance_sumsq,
        valid_maps=valid_maps,
        minimum_finite_stars=minimum_finite_stars,
        maximum_finite_stars=maximum_finite_stars,
        retained_stars=np.asarray(len(rows), dtype=np.int64),
        retained_maps=np.asarray(processed_maps, dtype=np.int64),
    )
    temporary.replace(path)
    return path


