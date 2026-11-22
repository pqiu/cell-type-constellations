"""Materialize the final 2,000-gene by 342,559-star PBMC matrix from source star maps."""
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from scripts.recompute_hierarchy_stars_h5ad import parent_design
from scripts.select_hierarchy_all_genes import ROOT, SOURCE, reconstruct_sample
from scripts.star_matrix_core import symbol_groups

BASE = ROOT / "results/hierarchy_global_reanalysis/final_132/pbmc_66_within_map_gene_selection_v1"
COHORT = ROOT / "results/hierarchy_global_reanalysis/final_132/pbmc_66_cohort_v1/star_columns.csv"
OUTPUT = BASE / "star_matrix_v1"
MATRIX = OUTPUT / "star_expression_2000_by_342559.float32.npy"
PROGRESS = OUTPUT / "progress.json"


def digest(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def atomic_json(path: Path, value: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def collapse_sample(dataset_id, sample_id, genes, groups, source_labels, incidence,
                    direct, small, old_h5, small_h5):
    present_positions = np.asarray([i for i, gene in enumerate(genes) if gene in groups], dtype=np.int64)
    present_genes = [genes[i] for i in present_positions]
    expanded_rows = np.concatenate([groups[g] for g in present_genes])
    expanded_names = [f"__native_row_{i}" for i in range(len(expanded_rows))]
    expanded_axis = dict(zip(expanded_names, map(int, expanded_rows)))
    counts = np.asarray([len(groups[g]) for g in present_genes], dtype=np.int64)
    starts = np.r_[0, np.cumsum(counts)[:-1]].astype(np.int64)
    means, parent_counts = reconstruct_sample(
        dataset_id, sample_id, expanded_names, expanded_axis, source_labels,
        incidence, direct, small, old_h5, small_h5,
    )
    finite = np.isfinite(means)
    totals = np.add.reduceat(np.where(finite, means, 0.0), starts, axis=1)
    observed = np.add.reduceat(finite.astype(np.int16), starts, axis=1)
    collapsed = np.divide(totals, observed, out=np.full(totals.shape, np.nan), where=observed > 0)
    return collapsed, parent_counts, present_positions


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    meta = pd.read_csv(COHORT)
    meta["matrix_column"] = np.arange(len(meta), dtype=np.int64)
    genes = (BASE / "selected_genes.txt").read_text(encoding="utf-8").splitlines()
    if len(meta) != 342_559 or meta.parent_id.nunique() != 66 or len(genes) != 2_000:
        raise ValueError("Final cohort or gene panel changed")
    if genes != sorted(set(genes)):
        raise ValueError("Selected genes must be sorted and unique")
    meta.to_csv(OUTPUT / "star_columns.csv", index=False)

    if PROGRESS.exists():
        progress = json.loads(PROGRESS.read_text(encoding="utf-8"))
        matrix = np.lib.format.open_memmap(MATRIX, mode="r+")
    else:
        matrix = np.lib.format.open_memmap(
            MATRIX, mode="w+", dtype=np.float32, shape=(len(genes), len(meta)), fortran_order=True
        )
        matrix[:] = np.nan
        matrix.flush()
        progress = {"completed_sources": [], "created_utc": datetime.now(timezone.utc).isoformat()}
        atomic_json(PROGRESS, progress)

    completed = set(progress["completed_sources"])
    source_labels, parents, _, incidence = parent_design(all_ancestors=True)
    parent_position = {name: i for i, name in enumerate(parents)}
    grouped = list(meta.groupby("dataset_id", sort=True))
    for number, (dataset_id, rows) in enumerate(grouped, 1):
        if dataset_id in completed:
            continue
        groups = symbol_groups(ROOT / "results/star_maps" / dataset_id / "genes.csv")
        curated = ROOT / "results/star_maps" / dataset_id
        hierarchy = SOURCE / dataset_id
        direct = pd.read_csv(curated / "stars.csv", usecols=["sample_id", "label", "n_cells", "status"])
        direct["direct_row"] = np.arange(len(direct), dtype=np.int64)
        small = pd.read_csv(hierarchy / "under20_direct_rows.csv")
        direct_path = curated / "expression.h5"
        if not direct_path.exists():
            direct_path = hierarchy / "native_direct_mean_basis.h5"
        with h5py.File(direct_path) as old_h5, h5py.File(hierarchy / "under20_direct_sums_native.h5") as small_h5:
            for sample_id, sample_rows in rows.groupby("sample_id", sort=False):
                collapsed, parent_counts, present = collapse_sample(
                    dataset_id, sample_id, genes, groups, source_labels, incidence,
                    direct, small, old_h5, small_h5,
                )
                positions = np.asarray([parent_position[x] for x in sample_rows.parent_id], dtype=np.int64)
                if not np.array_equal(parent_counts[positions], sample_rows.n_cells.to_numpy(dtype=np.int64)):
                    raise ValueError(f"Parent support changed: {dataset_id} {sample_id}")
                columns = sample_rows.matrix_column.to_numpy(dtype=np.int64)
                matrix[present[:, None], columns[None, :]] = collapsed[positions].astype(np.float32).T
        matrix.flush()
        progress["completed_sources"].append(dataset_id)
        atomic_json(PROGRESS, progress)
        if number % 25 == 0 or number == len(grouped):
            print(f"MATRIX {number}/{len(grouped)} sources", flush=True)

    means = np.empty(len(genes), dtype=np.float64)
    imputed = np.empty(len(genes), dtype=np.int64)
    for row in range(len(genes)):
        values = np.asarray(matrix[row], dtype=np.float64)
        finite = np.isfinite(values)
        if not finite.any():
            raise ValueError(f"Gene {genes[row]} has no observed values")
        means[row] = values[finite].mean()
        imputed[row] = int((~finite).sum())
        matrix[row, ~finite] = means[row]
    matrix.flush()

    stats = pd.read_csv(BASE / "gene_statistics.csv").set_index("gene").loc[genes].reset_index()
    gene_rows = pd.DataFrame({
        "matrix_row": np.arange(len(genes), dtype=np.int32), "gene": genes,
        "imputation_mean": means, "imputed_columns": imputed,
    }).merge(stats, on="gene", how="left", validate="one_to_one")
    gene_rows.to_csv(OUTPUT / "gene_rows.csv", index=False)
    if not np.isfinite(matrix).all():
        raise ValueError("Materialized matrix contains non-finite values")
    receipt = {
        "complete": True, "verified": True,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "matrix_shape": [2_000, len(meta)], "matrix_dtype": "float32",
        "genes": 2_000, "stars": len(meta), "maps": int(meta.map_index.nunique()),
        "types": int(meta.parent_id.nunique()), "sources": int(meta.dataset_id.nunique()),
        "all_values_finite": True, "matrix_sha256": digest(MATRIX),
        "star_columns_sha256": digest(OUTPUT / "star_columns.csv"),
        "gene_rows_sha256": digest(OUTPUT / "gene_rows.csv"),
    }
    atomic_json(OUTPUT / "completion.json", receipt)
    print(json.dumps(receipt, indent=2), flush=True)


if __name__ == "__main__":
    main()
