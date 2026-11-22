"""Build Figure 2 from the finalized 66-star raw and aligned matrices."""
from __future__ import annotations

import concurrent.futures
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist, pdist
import torch


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "results/hierarchy_global_reanalysis/final_132/pbmc_66_within_map_gene_selection_v1"
RAW_DIR = BASE / "star_matrix_v1"
ALIGN_DIR = BASE / "whole_constellation_alignment_v1"
REF_DIR = BASE / "reference_map_v1"
OUT = BASE / "figure2_whole_constellation_v1"
RAW_PATH = RAW_DIR / "star_expression_2000_by_342559.float32.npy"
ALIGNED_PATH = ALIGN_DIR / "aligned_expression_2000_by_342559.float32.npy"


def sha256(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def atomic_json(path: Path, value: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def pairwise_mean(values: torch.Tensor, block: int = 384) -> float:
    n = len(values)
    total = torch.zeros((), dtype=torch.float64, device=values.device)
    for start in range(0, n, block):
        stop = min(start + block, n)
        distances = torch.cdist(values[start:stop], values)
        distances[
            torch.arange(stop - start, device=values.device),
            torch.arange(start, stop, device=values.device),
        ] = 0
        total += distances.sum(dtype=torch.float64)
    return float((total / (n * (n - 1))).item())


def score_worker(device_index: int, assignments: list[dict]) -> list[dict]:
    torch.cuda.set_device(device_index)
    raw = np.load(RAW_PATH, mmap_mode="r")
    aligned = np.load(ALIGNED_PATH, mmap_mode="r")
    records = []
    for position, assignment in enumerate(assignments, 1):
        columns = np.asarray(assignment["columns"], dtype=np.int64)
        raw_tensor = torch.as_tensor(
            np.ascontiguousarray(raw[:, columns].T), dtype=torch.float32,
            device=f"cuda:{device_index}",
        )
        aligned_tensor = torch.as_tensor(
            np.ascontiguousarray(aligned[:, columns].T), dtype=torch.float32,
            device=f"cuda:{device_index}",
        )
        raw_mean = pairwise_mean(raw_tensor)
        aligned_mean = pairwise_mean(aligned_tensor)
        rng = np.random.default_rng(20261114 + int(assignment["reference_column"]))
        subset = rng.choice(len(columns), min(24, len(columns)), replace=False)
        direct = pdist(np.asarray(raw[:, columns[subset]].T, dtype=np.float64))
        gpu = torch.pdist(raw_tensor[torch.as_tensor(subset, device=raw_tensor.device)]).double().cpu().numpy()
        records.append({
            "reference_column": assignment["reference_column"],
            "parent_id": assignment["parent_id"],
            "cell_type_name": assignment["cell_type_name"],
            "maps": len(columns),
            "map_pairs": len(columns) * (len(columns) - 1) // 2,
            "raw_distance": raw_mean,
            "aligned_distance": aligned_mean,
            "gpu_subset_max_distance_error": float(np.max(np.abs(direct - gpu))) if len(direct) else 0.0,
            "gpu_device": device_index,
        })
        print(
            f"GPU {device_index}: {position}/{len(assignments)} {assignment['cell_type_name']} "
            f"raw={raw_mean:.5f} aligned={aligned_mean:.5f}", flush=True,
        )
        del raw_tensor, aligned_tensor
        torch.cuda.empty_cache()
    return records


def render_figure(table: pd.DataFrame) -> None:
    """Render the verified Figure 2 table without recomputing distances."""
    maximum = float(max(
        table.within_map_distance.max(), table.raw_distance.max(), table.aligned_distance.max()
    ) * 1.07)
    plt.rcParams.update({
        "font.size": 11, "axes.spines.top": False, "axes.spines.right": False,
        "svg.fonttype": "none", "font.family": "DejaVu Sans",
    })
    fig, axes = plt.subplots(1, 2, figsize=(13.8, 7.2), sharex=True, sharey=True)
    for ax, column, title, letter in zip(
        axes, ["raw_distance", "aligned_distance"],
        ["Before alignment", "After one whole-constellation shift"], ["A", "B"],
    ):
        ax.plot([0, maximum], [0, maximum], color="#999999", lw=1, ls="--", zorder=0)
        ax.scatter(table.within_map_distance, table[column], s=31, color="#256799",
                   alpha=0.76, edgecolors="none", zorder=2)
        ax.set(xlim=(0, maximum), ylim=(0, maximum)); ax.set_aspect("equal", adjustable="box")
        ax.set_title(title, loc="left", pad=12, fontsize=14)
        ax.text(maximum * 0.88, maximum * 0.93, "y = x", ha="center", va="center",
                fontsize=10, color="#777777", rotation=45)
        ax.text(-0.09, 1.035, letter, transform=ax.transAxes, fontsize=15, fontweight="bold")
    axes[1].tick_params(labelleft=True)
    fig.suptitle("Entire 66-star map as one constellation", fontsize=17, y=0.975)
    fig.supxlabel("Mean distance to other types within the same map", y=0.030, fontsize=12)
    fig.supylabel("Mean distance of the same type across maps", x=0.045, fontsize=12)
    fig.subplots_adjust(left=0.085, right=0.985, bottom=0.13, top=0.87, wspace=0.16)
    for extension in ["png", "svg", "pdf"]:
        fig.savefig(OUT / f"figure2_whole_constellation.{extension}", dpi=240)
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    raw = np.load(RAW_PATH, mmap_mode="r")
    aligned = np.load(ALIGNED_PATH, mmap_mode="r")
    stars = pd.read_csv(RAW_DIR / "star_columns.csv")
    types = pd.read_csv(REF_DIR / "reference_star_columns.csv")
    if raw.shape != aligned.shape or raw.shape != (2_000, 342_559) or len(types) != 66:
        raise ValueError("Unexpected Figure 2 input inventory")

    lookup = dict(zip(types.parent_id, types.reference_column))
    sums = np.zeros(len(types), dtype=np.float64)
    counts = np.zeros(len(types), dtype=np.int64)
    for completed, (_, group) in enumerate(stars.groupby("map_index", sort=True), 1):
        columns = group.matrix_column.to_numpy(dtype=np.int64)
        values = np.asarray(raw[:, columns].T, dtype=np.float64)
        distances = cdist(values, values, metric="euclidean")
        means = distances.sum(axis=1) / (len(columns) - 1)
        type_columns = np.fromiter(
            (lookup[x] for x in group.parent_id), dtype=np.int64, count=len(group)
        )
        np.add.at(sums, type_columns, means)
        np.add.at(counts, type_columns, 1)
        if completed % 2_000 == 0:
            print(f"WITHIN {completed:,}/14,230 maps", flush=True)
    within = types.copy()
    within["within_map_distance"] = sums / counts
    within["within_map_count"] = counts
    within.to_csv(OUT / "within_map_per_type.csv", index=False)

    assignments = []
    for row in types.itertuples(index=False):
        columns = stars.loc[stars.parent_id.eq(row.parent_id), "matrix_column"].tolist()
        assignments.append({
            "reference_column": int(row.reference_column), "parent_id": row.parent_id,
            "cell_type_name": row.cell_type_name, "columns": columns,
            "cost": len(columns) ** 2,
        })
    devices = [[], []]
    loads = [0, 0]
    for assignment in sorted(assignments, key=lambda x: x["cost"], reverse=True):
        device = int(np.argmin(loads))
        devices[device].append(assignment)
        loads[device] += assignment["cost"]
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        futures = [pool.submit(score_worker, device, devices[device]) for device in range(2)]
        records = []
        for future in futures:
            records.extend(future.result())
    across = pd.DataFrame(records).sort_values("reference_column")
    table = types.merge(
        within[["parent_id", "within_map_distance", "within_map_count"]], on="parent_id"
    ).merge(
        across.drop(columns=["cell_type_name", "reference_column"]), on="parent_id"
    )
    table["reduction_percent"] = 100 * (1 - table.aligned_distance / table.raw_distance)
    table["raw_over_within"] = table.raw_distance / table.within_map_distance
    table["aligned_over_within"] = table.aligned_distance / table.within_map_distance
    table.to_csv(OUT / "per_type_figure2.csv", index=False, float_format="%.10g")

    render_figure(table)

    summary = {
        "complete": True, "verified": True,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "maps": int(stars.map_index.nunique()), "stars": len(stars),
        "types": len(table), "genes": raw.shape[0],
        "mean_reduction_percent": float(table.reduction_percent.mean()),
        "median_reduction_percent": float(table.reduction_percent.median()),
        "types_improved": int((table.reduction_percent > 0).sum()),
        "types_below_diagonal_before": int((table.raw_distance < table.within_map_distance).sum()),
        "types_below_diagonal_after": int((table.aligned_distance < table.within_map_distance).sum()),
        "median_raw_over_within": float(table.raw_over_within.median()),
        "median_aligned_over_within": float(table.aligned_over_within.median()),
        "maximum_gpu_subset_distance_error": float(table.gpu_subset_max_distance_error.max()),
        "raw_matrix_sha256": sha256(RAW_PATH), "aligned_matrix_sha256": sha256(ALIGNED_PATH),
        "per_type_sha256": sha256(OUT / "per_type_figure2.csv"),
        "figure_png_sha256": sha256(OUT / "figure2_whole_constellation.png"),
        "limitation": "Pooled descriptive map-view analysis; repeated maps, donors and studies are not independent and the inclusive shift partly uses each evaluated star.",
    }
    atomic_json(OUT / "completion.json", summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    if "--plot-only" in sys.argv:
        render_figure(pd.read_csv(OUT / "per_type_figure2.csv"))
        completion_path = OUT / "completion.json"
        summary = json.loads(completion_path.read_text(encoding="utf-8"))
        summary["figure_rendered_utc"] = datetime.now(timezone.utc).isoformat()
        summary["figure_png_sha256"] = sha256(OUT / "figure2_whole_constellation.png")
        atomic_json(completion_path, summary)
    else:
        main()
