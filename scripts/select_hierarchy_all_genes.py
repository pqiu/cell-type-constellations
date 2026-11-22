"""Select training-only genes from the validated, screened hierarchy stars.

Native profiles are reconstructed from the saved direct-label means and
under-20-cell sums. This does not reopen counts or alter the unfiltered stars.
"""
from __future__ import annotations

import json
import hashlib
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from scripts.inspect_pooled_star_axes import ROOT, eligible_outputs, symbol_column
from scripts.recompute_hierarchy_stars_h5ad import parent_design
from src.hierarchy_stars import pooled_parent_means

SOURCE = ROOT / "results/hierarchy_all_star_maps"
OUTPUT = ROOT / "results/hierarchy_global_reanalysis/final_132"
SEED = 41
SAMPLED_STARS = 20_000
GENES = 2_000


def unique_symbol_axis(path: Path) -> dict[str, int]:
    genes = pd.read_csv(path, keep_default_na=False)
    field = symbol_column(genes.columns)
    names = genes[field].astype(str).str.strip()
    counts = names.value_counts()
    return {name: i for i, name in enumerate(names) if name and counts[name] == 1}


def screened_inventory() -> tuple[pd.DataFrame, list[str], dict[str, dict[str, int]]]:
    receipt = json.loads((SOURCE / "full_validation.json").read_text())
    if receipt["sources"] != 961 or receipt["sample_views"] != 16_524:
        raise ValueError("Full hierarchy validation is missing")
    coverage = pd.read_csv(SOURCE / "all_ancestor_coverage.csv")
    retained = set(coverage.loc[coverage.retained_after_20_all_maps, "aggregate_id"])
    if len(retained) != 132:
        raise ValueError("Across-collection type screen changed")
    split_rows = pd.read_csv(ROOT / "results/pooled_star_translation/maps.csv",
                             usecols=["study_id", "split"]).drop_duplicates("study_id")
    split = split_rows.set_index("study_id").split.to_dict()
    axes: dict[str, dict[str, int]] = {}
    core: set[str] | None = None
    blocks = []
    for dataset_id in eligible_outputs():
        curated = ROOT / "results/star_maps" / dataset_id
        run = json.loads((curated / "run.json").read_text())
        study = run["configuration"]["study_id"]
        if study not in split:
            raise ValueError(f"Study lacks frozen split: {study}")
        axis = unique_symbol_axis(curated / "genes.csv")
        axes[dataset_id] = axis
        core = set(axis) if core is None else core & set(axis)
        stars = pd.read_csv(SOURCE / dataset_id / "stars.csv",
                            usecols=["sample_id", "donor_uid", "parent_id", "n_cells",
                                     "status", "expression_row"])
        keep = stars.parent_id.isin(retained) & stars.n_cells.ge(20)
        block = stars.loc[keep].copy()
        if not block.status.eq("supported").all():
            raise ValueError(f"Support screen differs: {dataset_id}")
        block["dataset_id"] = dataset_id
        block["study_id"] = study
        block["split"] = split[study]
        blocks.append(block)
    meta = pd.concat(blocks, ignore_index=True)
    if len(meta) != 875_669 or core is None or len(core) != 6_870:
        raise ValueError(f"Filtered inventory or shared axis changed: {len(meta)}, {len(core or [])}")
    return meta, sorted(core), axes


def reconstruct_sample(dataset_id: str, sample_id: str,
                       genes: list[str], axis: dict[str, int], source_labels: list[str],
                       incidence, direct: pd.DataFrame, small: pd.DataFrame,
                       old_h5: h5py.File, small_h5: h5py.File) -> tuple[np.ndarray, np.ndarray]:
    direct_sample = direct.loc[direct.sample_id.eq(sample_id)].set_index("label").reindex(source_labels)
    if direct_sample.n_cells.isna().any():
        raise ValueError(f"Incomplete direct-label counts: {dataset_id} {sample_id}")
    counts = direct_sample.n_cells.to_numpy(dtype=np.int64)
    columns = np.asarray([axis[g] for g in genes], dtype=np.int64)
    sums = np.zeros((len(source_labels), len(genes)), dtype=np.float64)
    supported = np.flatnonzero(direct_sample.status.eq("supported").to_numpy())
    if len(supported):
        rows = direct_sample.iloc[supported].direct_row.to_numpy(dtype=np.int64)
        if not np.all(np.diff(rows) > 0):
            raise ValueError("Direct-label mean rows are not ordered")
        means = old_h5["norm10k_log1p_mean"][rows, :][:, columns]
        sums[supported] = means.astype(np.float64) * counts[supported, None]
    low = direct_sample.status.eq("insufficient_support").to_numpy()
    if low.any():
        low_labels = np.asarray(source_labels)[low]
        small_sample = small.loc[small.sample_id.eq(sample_id)].set_index("label")
        rows = small_sample.loc[low_labels, "basis_row"].to_numpy(dtype=np.int64)
        if not np.all(np.diff(rows) > 0):
            raise ValueError("Under-20 basis rows are not ordered")
        sums[np.flatnonzero(low)] = small_h5["normalized_log_sum"][rows, :][:, columns]
    # Reconstruct only requested parents; the incidence matrix has the same
    # sorted direct-label and ancestor ordering as the saved unfiltered run.
    means, parent_counts = pooled_parent_means(sums, counts, incidence, min_cells=1)
    return means, parent_counts


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    meta, core, axes = screened_inventory()
    source_labels, all_parents, _, incidence = parent_design(all_ancestors=True)
    parent_pos = {name: i for i, name in enumerate(all_parents)}
    train = meta.loc[meta.split.eq("train")].copy()
    # Repeated views of a donor/type do not get multiple chances to determine
    # the gene panel. The deterministic first view is its representative.
    representatives = train.drop_duplicates(["study_id", "donor_uid", "parent_id"])
    rng = np.random.default_rng(SEED)
    chosen = np.sort(rng.choice(len(representatives),
                                min(SAMPLED_STARS, len(representatives)), replace=False))
    sampled = representatives.iloc[chosen].reset_index(drop=True)
    if len(sampled) != SAMPLED_STARS:
        raise ValueError("Too few training donor/type stars")
    sampled.to_csv(OUTPUT / "gene_selection_stars.csv", index=False)
    native = np.memmap(OUTPUT / "gene_selection_native.float32", mode="w+",
                       dtype="float32", shape=(len(sampled), len(core)))
    for k, (dataset_id, group) in enumerate(sampled.groupby("dataset_id", sort=False), 1):
        curated = ROOT / "results/star_maps" / dataset_id
        hierarchy = SOURCE / dataset_id
        direct = pd.read_csv(curated / "stars.csv",
                             usecols=["sample_id", "label", "n_cells", "status"])
        direct["direct_row"] = np.arange(len(direct), dtype=np.int64)
        small = pd.read_csv(hierarchy / "under20_direct_rows.csv")
        old_path = curated / "expression.h5"
        if not old_path.exists():
            old_path = hierarchy / "native_direct_mean_basis.h5"
            receipt = json.loads((hierarchy / "retired_direct_label_expression.json").read_text())
            with old_path.open("rb") as handle:
                basis_hash = hashlib.file_digest(handle, "sha256").hexdigest()
            if receipt["basis_sha256"] != basis_hash:
                raise ValueError(f"Retired direct basis hash changed: {dataset_id}")
        with h5py.File(old_path) as old_h5, h5py.File(hierarchy / "under20_direct_sums_native.h5") as small_h5:
            for sample_id, rows in group.groupby("sample_id", sort=False):
                means, counts = reconstruct_sample(dataset_id, sample_id,
                                                   core, axes[dataset_id], source_labels,
                                                   incidence, direct, small, old_h5, small_h5)
                for row in rows.itertuples():
                    j = parent_pos[row.parent_id]
                    if counts[j] != row.n_cells:
                        raise ValueError(f"Reconstructed parent support changed: {dataset_id} {sample_id} {row.parent_id}")
                    native[row.Index] = means[j].astype(np.float32)
        if k % 50 == 0:
            print(f"RECONSTRUCTED {k} sources for training-only gene selection", flush=True)
    native.flush()
    finite = np.isfinite(native).all(axis=0)
    if finite.sum() < GENES:
        raise ValueError(f"Only {finite.sum()} complete genes in training sample")
    variance = np.var(native[:, finite], axis=0, dtype=np.float64)
    available = np.asarray(core)[finite]
    ranking = np.lexsort((available, -variance))
    selected = sorted(available[ranking[:GENES]].tolist())
    (OUTPUT / "selected_genes.txt").write_text("\n".join(selected) + "\n")
    report = {"sources": 961, "filtered_stars": len(meta), "retained_types": 132,
              "shared_unique_genes": len(core), "training_donor_type_stars": len(representatives),
              "sampled_training_stars": len(sampled), "seed": SEED,
              "complete_sampled_genes": int(finite.sum()), "selected_genes": len(selected),
              "split_source": "results/pooled_star_translation/maps.csv, frozen study split",
              "selection": "Top variance of 20,000 seeded training-study donor/type representative native hierarchy centroids",
              "training_only": True}
    (OUTPUT / "gene_selection_report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
