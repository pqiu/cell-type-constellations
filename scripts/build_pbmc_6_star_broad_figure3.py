"""Build Figure 3 with six existing broad PBMC parent stars."""
from __future__ import annotations

import hashlib
import importlib.metadata
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

os.environ["NUMBA_NUM_THREADS"] = "8"

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from threadpoolctl import threadpool_limits
import umap
from pynndescent import rp_trees


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "results/hierarchy_global_reanalysis/final_132/pbmc_71_within_map_gene_selection_v3"
RAW_DIR = BASE / "star_matrix_v1"
REF_DIR = BASE / "reference_map_v1"
OUT = BASE / "figure3_six_star_broad_v1"
RAW_MATRIX_NAME = "star_expression_2000_by_345670.float32.npy"
REFERENCE_MATRIX_NAME = "reference_expression_2000_by_71.float64.npy"
SEED = 41
MIN_MEMBERS = 3
BROAD = {
    "CL:0000236": "B cell",
    "CL:0000624": "CD4 T cell",
    "CL:0000625": "CD8 T cell",
    "CL:0000623": "Natural killer cell",
    "CL:0000576": "Monocyte",
    "CL:0000451": "Dendritic cell",
}
COLORS = ["#1f77b4", "#74c476", "#006d2c", "#9467bd", "#8c564b", "#202020"]


def sha256(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def write_json(path: Path, value: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def progress(name: str, start: float, **extra) -> None:
    write_json(OUT / "progress.json", {"stage": name, "elapsed_seconds": time.monotonic() - start, **extra})
    print(name, flush=True)


def render_figure(matched: pd.DataFrame, coordinates: dict[str, np.ndarray]) -> None:
    """Render Figure 3 from saved, verified coordinates."""
    color_map = dict(zip(BROAD, COLORS))
    labels = matched.parent_id.to_numpy()
    point_colors = np.asarray([color_map[x] for x in labels])
    order = np.random.default_rng(SEED).permutation(len(matched))
    handles = [
        plt.Line2D([0], [0], marker="o", linestyle="", color=color_map[x], label=BROAD[x], markersize=7)
        for x in BROAD
    ]
    plt.rcParams.update({"font.family": "DejaVu Sans", "svg.fonttype": "none"})
    fig, axes = plt.subplots(1, 2, figsize=(15, 6.6))
    for ax, mode, letter, title in zip(
        axes, ["raw", "aligned"], ["A", "B"],
        ["Before alignment", "After 6-star constellation alignment"],
    ):
        xy = coordinates[mode]
        ax.scatter(xy[order, 0], xy[order, 1], s=2, c=point_colors[order], alpha=0.45,
                   linewidths=0, rasterized=True)
        ax.set_title(title, loc="left", fontsize=14, pad=10)
        ax.text(-0.08, 1.025, letter, transform=ax.transAxes, fontsize=15, fontweight="bold")
        ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
        ax.set_xticks([]); ax.set_yticks([]); ax.spines[["top", "right"]].set_visible(False)
        ax.set_box_aspect(1)
    fig.suptitle("A broad 6-star constellation spanning 5 immune lineages", fontsize=17, y=0.975)
    fig.legend(handles=handles, loc="center left", ncol=1, bbox_to_anchor=(0.835, 0.49), frameon=False, fontsize=11)
    fig.subplots_adjust(left=0.055, right=0.82, top=0.86, bottom=0.11, wspace=0.16)
    for extension in ["png", "svg"]:
        fig.savefig(OUT / f"figure3_six_star_broad.{extension}", dpi=240)
    plt.close(fig)


def main() -> None:
    start = time.monotonic()
    OUT.mkdir(parents=True, exist_ok=True)
    progress("select_maps", start)
    raw_matrix = np.load(RAW_DIR / RAW_MATRIX_NAME, mmap_mode="r")
    reference = np.load(REF_DIR / REFERENCE_MATRIX_NAME, mmap_mode="r")
    stars = pd.read_csv(RAW_DIR / "star_columns.csv")
    reference_stars = pd.read_csv(REF_DIR / "reference_star_columns.csv")
    reference_lookup = dict(zip(reference_stars.parent_id, reference_stars.reference_column))
    if not set(BROAD).issubset(reference_lookup):
        raise ValueError(f"Missing broad parents: {sorted(set(BROAD) - set(reference_lookup))}")

    panel = stars.loc[stars.parent_id.isin(BROAD)].copy()
    counts = panel.groupby("map_index").size()
    qualifying_maps = counts.index[counts >= MIN_MEMBERS].to_numpy(dtype=np.int64)
    matched = panel.loc[panel.map_index.isin(qualifying_maps)].copy()
    matched.sort_values("matrix_column", kind="stable", inplace=True)
    matched.reset_index(drop=True, inplace=True)
    matched.insert(0, "matched_profile", np.arange(len(matched)))
    matched["display_name"] = matched.parent_id.map(BROAD)
    matched.to_csv(OUT / "matched_stars.csv", index=False)
    pd.DataFrame({
        "map_index": qualifying_maps,
        "observed_panel_members": counts.loc[qualifying_maps].to_numpy(dtype=np.int16),
    }).to_csv(OUT / "qualifying_maps.csv", index=False)
    members = reference_stars.loc[reference_stars.parent_id.isin(BROAD)].copy()
    members["display_name"] = members.parent_id.map(BROAD)
    members.to_csv(OUT / "members.csv", index=False)

    columns = matched.matrix_column.to_numpy(dtype=np.int64)
    raw = np.ascontiguousarray(raw_matrix[:, columns].T)
    if raw.shape != (len(matched), 2_000) or not np.isfinite(raw).all():
        raise ValueError("Invalid broad-panel raw matrix")

    progress("align_one_shift_per_map", start, maps=len(qualifying_maps), stars=len(matched))
    aligned = np.empty_like(raw)
    maximum_constraint_error = 0.0
    for completed, (map_index, positions) in enumerate(matched.groupby("map_index", sort=True).indices.items(), 1):
        positions = np.asarray(positions, dtype=np.int64)
        ref_columns = np.fromiter(
            (reference_lookup[x] for x in matched.iloc[positions].parent_id),
            dtype=np.int64,
            count=len(positions),
        )
        shift = np.mean(raw[positions].astype(np.float64) - reference[:, ref_columns].T, axis=0)
        aligned[positions] = raw[positions] - shift
        error = np.max(np.abs(np.mean(aligned[positions].astype(np.float64) - reference[:, ref_columns].T, axis=0)))
        maximum_constraint_error = max(maximum_constraint_error, float(error))
        if completed % 2_000 == 0:
            print(f"ALIGN {completed:,}/{len(qualifying_maps):,}", flush=True)
    if maximum_constraint_error > 5e-7 or not np.isfinite(aligned).all():
        raise ValueError(f"Broad alignment failed: {maximum_constraint_error}")
    np.save(OUT / "aligned_profiles.float32.npy", aligned)

    progress("fit_shared_raw_pca", start)
    with threadpool_limits(limits=8):
        pca = PCA(n_components=50, svd_solver="randomized", random_state=SEED)
        raw_scores = pca.fit_transform(raw).astype(np.float32)
        aligned_scores = pca.transform(aligned).astype(np.float32)
    np.save(OUT / "raw_pca_scores.npy", raw_scores)
    np.save(OUT / "aligned_scores_in_raw_pca.npy", aligned_scores)
    np.savez(
        OUT / "shared_raw_fitted_pca.npz",
        components=pca.components_, mean=pca.mean_,
        explained_variance=pca.explained_variance_,
        explained_variance_ratio=pca.explained_variance_ratio_,
    )

    rp_trees.rptree_leaf_array_parallel = lambda forest: [
        rp_trees.get_leaves_from_tree(tree, max(tree.leaf_size for tree in forest)) for tree in forest
    ]
    coordinates = {}
    for mode, values in [("raw", raw_scores), ("aligned", aligned_scores)]:
        progress(f"fit_{mode}_umap", start)
        reducer = umap.UMAP(
            n_neighbors=40, min_dist=0.2, n_components=2, metric="euclidean",
            random_state=SEED, transform_seed=SEED, n_jobs=1,
        )
        with threadpool_limits(limits=8):
            xy = reducer.fit_transform(values).astype(np.float32)
        if xy.shape != (len(matched), 2) or not np.isfinite(xy).all():
            raise ValueError(f"Invalid {mode} UMAP")
        np.save(OUT / f"{mode}_umap_coordinates.npy", xy)
        coordinates[mode] = xy
        matched[f"{mode}_umap1"] = xy[:, 0]
        matched[f"{mode}_umap2"] = xy[:, 1]
    matched.to_csv(OUT / "umap_coordinates.csv", index=False)

    progress("render", start)
    render_figure(matched, coordinates)

    report = {
        "complete": True,
        "verified": True,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "panel_members": BROAD,
        "members_are_existing_broad_parent_stars": True,
        "minimum_observed_members_per_map": MIN_MEMBERS,
        "qualifying_maps": len(qualifying_maps),
        "studies": int(matched.study_id.nunique()),
        "stars_per_panel": len(matched),
        "identical_stars_in_both_panels": True,
        "alignment": "one equal-star shift per map from every observed broad-panel member; no target-aware exclusions",
        "maximum_aligned_map_mean_displacement_error": maximum_constraint_error,
        "shared_pca": "50 components fitted to unaligned matched stars only; same basis transforms aligned stars",
        "raw_pca_explained_variance": float(pca.explained_variance_ratio_.sum()),
        "umap": {"neighbors": 40, "minimum_distance": 0.2, "metric": "euclidean", "seed": SEED, "independent_fits": True},
        "gene_panel_sha256": sha256(BASE / "selected_genes.txt"),
        "figure_png_sha256": sha256(OUT / "figure3_six_star_broad.png"),
        "raw_coordinates_sha256": sha256(OUT / "raw_umap_coordinates.npy"),
        "aligned_coordinates_sha256": sha256(OUT / "aligned_umap_coordinates.npy"),
        "versions": {package: importlib.metadata.version(package) for package in [
            "numpy", "pandas", "scikit-learn", "umap-learn", "matplotlib"
        ]},
        "elapsed_seconds": time.monotonic() - start,
    }
    write_json(OUT / "completion.json", report)
    progress("complete", start, maps=len(qualifying_maps), stars=len(matched))
    print(json.dumps(report, indent=2), flush=True)


def plot_only() -> None:
    saved_matched = pd.read_csv(OUT / "matched_stars.csv")
    saved_coordinates = {
        "raw": np.load(OUT / "raw_umap_coordinates.npy"),
        "aligned": np.load(OUT / "aligned_umap_coordinates.npy"),
    }
    render_figure(saved_matched, saved_coordinates)
    completion_path = OUT / "completion.json"
    report = json.loads(completion_path.read_text(encoding="utf-8"))
    report["figure_rendered_utc"] = datetime.now(timezone.utc).isoformat()
    report["figure_png_sha256"] = sha256(OUT / "figure3_six_star_broad.png")
    write_json(completion_path, report)


if __name__ == "__main__":
    if "--plot-only" in sys.argv:
        plot_only()
    else:
        main()
