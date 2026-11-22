"""Select 2,000 genes by summed within-map variance in the frozen 66-star cohort."""
from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
sys.path.insert(0, str(SCRIPT_DIR.parent))

from scripts.recompute_hierarchy_stars_h5ad import parent_design
import scripts.within_map_variance_core as legacy
from scripts.select_hierarchy_all_genes import ROOT


FINAL = ROOT / "results/hierarchy_global_reanalysis/final_132"
COHORT = FINAL / "pbmc_66_cohort_v1"
OUTPUT = FINAL / "pbmc_66_within_map_gene_selection_v1"
PARTIALS = OUTPUT / "source_partials"
SELECTED_GENES = 2_000


def atomic_json(path: Path, value: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def main() -> None:
    start = time.monotonic()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    meta = pd.read_csv(COHORT / "star_columns.csv")
    availability = pd.read_csv(COHORT / "gene_availability.csv")
    genes = sorted(availability.loc[availability.passes_90_percent, "gene"].tolist())
    if len(genes) != 17_900 or len(set(genes)) != len(genes):
        raise ValueError(f"Expected 17,900 eligible genes, found {len(genes)}")
    (OUTPUT / "eligible_genes.txt").write_text("\n".join(genes) + "\n", encoding="utf-8")

    legacy.PARTIALS = PARTIALS
    source_labels, parents, _, incidence = parent_design(all_ancestors=True)
    variance_sum = np.zeros(len(genes), dtype=np.float64)
    variance_sumsq = np.zeros(len(genes), dtype=np.float64)
    valid_maps = np.zeros(len(genes), dtype=np.int64)
    minimum_finite_stars = np.full(len(genes), np.iinfo(np.int16).max, dtype=np.int16)
    maximum_finite_stars = np.zeros(len(genes), dtype=np.int16)
    grouped = list(meta.groupby("dataset_id", sort=True))

    for number, (dataset_id, rows) in enumerate(grouped, 1):
        path = legacy.source_partial(
            dataset_id, rows, genes, source_labels, parents, incidence
        )
        with np.load(path, allow_pickle=False) as saved:
            indices = saved["gene_indices"]
            variance_sum[indices] += saved["variance_sum"]
            variance_sumsq[indices] += saved["variance_sumsq"]
            valid_maps[indices] += saved["valid_maps"]
            has_maps = saved["valid_maps"] > 0
            minimum_finite_stars[indices[has_maps]] = np.minimum(
                minimum_finite_stars[indices[has_maps]],
                saved["minimum_finite_stars"][has_maps],
            )
            maximum_finite_stars[indices] = np.maximum(
                maximum_finite_stars[indices], saved["maximum_finite_stars"]
            )
        if number % 10 == 0 or number == len(grouped):
            atomic_json(
                OUTPUT / "progress.json",
                {
                    "stage": "source_partials",
                    "completed_sources": number,
                    "total_sources": len(grouped),
                    "elapsed_seconds": time.monotonic() - start,
                },
            )
            print(f"WITHIN_MAP_VARIANCE {number}/{len(grouped)} sources", flush=True)

    if (valid_maps <= 0).any() or not np.isfinite(variance_sum).all():
        raise ValueError("A gene lacks valid map variances")
    mean_map_variance = variance_sum / valid_maps
    variance_of_map_variance = np.maximum(
        variance_sumsq / valid_maps - np.square(mean_map_variance), 0.0
    )
    order = np.lexsort((np.asarray(genes), -variance_sum))
    ranks = np.empty(len(genes), dtype=np.int64)
    ranks[order] = np.arange(1, len(genes) + 1)
    selected_ranked = [genes[i] for i in order[:SELECTED_GENES]]
    selected = set(selected_ranked)
    availability = availability.set_index("gene")

    statistics = pd.DataFrame(
        {
            "gene": genes,
            "available_columns": [availability.loc[g, "available_columns"] for g in genes],
            "available_fraction": [availability.loc[g, "available_fraction"] for g in genes],
            "maps_with_variance": valid_maps,
            "map_coverage_fraction": valid_maps / meta.map_index.nunique(),
            "sum_within_map_population_variance": variance_sum,
            "mean_within_map_population_variance": mean_map_variance,
            "sd_within_map_population_variance": np.sqrt(variance_of_map_variance),
            "minimum_finite_stars_in_scored_map": minimum_finite_stars,
            "maximum_finite_stars_in_scored_map": maximum_finite_stars,
            "within_map_variance_rank": ranks,
            "selected": [g in selected for g in genes],
        }
    ).sort_values(["within_map_variance_rank", "gene"])
    statistics.to_csv(OUTPUT / "gene_statistics.csv", index=False)
    (OUTPUT / "selected_genes_ranked.txt").write_text(
        "\n".join(selected_ranked) + "\n", encoding="utf-8"
    )
    (OUTPUT / "selected_genes.txt").write_text(
        "\n".join(sorted(selected_ranked)) + "\n", encoding="utf-8"
    )

    prior = set(
        (FINAL / "pbmc_71_within_map_gene_selection_v3/selected_genes.txt")
        .read_text(encoding="utf-8")
        .splitlines()
    )
    completion = {
        "complete": True,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "types": int(meta.parent_id.nunique()),
        "stars": int(len(meta)),
        "maps": int(meta.map_index.nunique()),
        "sources": int(meta.dataset_id.nunique()),
        "eligible_genes": len(genes),
        "selected_genes": SELECTED_GENES,
        "selected_overlap_with_prior_71_star_panel": len(selected & prior),
        "selected_new_relative_to_prior_71_star_panel": len(selected - prior),
        "elapsed_seconds": time.monotonic() - start,
    }
    atomic_json(OUTPUT / "completion.json", completion)
    print(json.dumps(completion, indent=2), flush=True)


if __name__ == "__main__":
    main()
