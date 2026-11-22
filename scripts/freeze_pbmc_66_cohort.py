"""Freeze the finalized 66-type PBMC cohort and its gene availability.

The cohort starts from the established 50-cell/50-map support filter, removes
the five manually reviewed anatomically incompatible labels, and reapplies the
minimum of 10 retained stars per sample-level map.
"""
from __future__ import annotations

from collections import Counter
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.inspect_pooled_star_axes import symbol_column
from scripts.select_hierarchy_all_genes import ROOT


FINAL = ROOT / "results/hierarchy_global_reanalysis/final_132"
FILTER = FINAL / "vocabulary_filter_50_cells_50_maps"
OUTPUT = FINAL / "pbmc_66_cohort_v1"
MIN_STARS_PER_MAP = 10
MIN_GENE_FRACTION = 0.90
EXCLUDED_TYPES = {
    "CL:1001603": "lung macrophage",
    "CL:0000844": "germinal center B cell",
    "CL:0000893": "thymocyte",
    "CL:0000809": "double-positive, alpha-beta thymocyte",
    "CL:0002489": "double negative thymocyte",
}


def source_symbols(dataset_id: str) -> set[str]:
    genes = pd.read_csv(
        ROOT / "results/star_maps" / dataset_id / "genes.csv",
        keep_default_na=False,
    )
    field = symbol_column(genes.columns)
    return set(genes[field].astype(str).str.strip()) - {""}


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    meta = pd.read_csv(FINAL / "star_rows.csv")
    supported_rows = np.load(FILTER / "supported_star_rows.npy")
    supported = meta.iloc[supported_rows].copy()
    if len(supported) != 358_928 or supported.parent_id.nunique() != 71:
        raise ValueError("The established 71-type support-filtered cohort changed")

    removed = supported.loc[supported.parent_id.isin(EXCLUDED_TYPES)].copy()
    retained = supported.loc[~supported.parent_id.isin(EXCLUDED_TYPES)].copy()
    pre_map_filter_stars = len(retained)
    pre_map_filter_maps = int(retained.map_index.nunique())

    map_counts = retained.groupby("map_index", sort=False).size()
    retained_map_indices = set(map_counts.loc[map_counts.ge(MIN_STARS_PER_MAP)].index)
    dropped_map_indices = set(map_counts.index) - retained_map_indices
    retained = retained.loc[retained.map_index.isin(retained_map_indices)].copy()
    retained.sort_values("star_row", inplace=True)
    retained.reset_index(drop=True, inplace=True)
    retained["matrix_column"] = np.arange(len(retained), dtype=np.int64)

    if retained.parent_id.nunique() != 66 or not retained.n_cells.ge(50).all():
        raise ValueError("Final 66-type cohort failed its support checks")
    if not retained.status.eq("supported").all():
        raise ValueError("Final cohort contains an unsupported star")

    type_names = pd.read_csv(FILTER / "type_eligibility.csv").set_index("parent_id")["name"]
    retained["cell_type_name"] = retained.parent_id.map(type_names)
    retained.to_csv(OUTPUT / "star_columns.csv", index=False)

    removal_counts = (
        removed.groupby("parent_id", sort=True)
        .size()
        .rename("stars_removed_before_map_filter")
        .reset_index()
    )
    removal_counts["cell_type_name"] = removal_counts.parent_id.map(EXCLUDED_TYPES)
    removal_counts.to_csv(OUTPUT / "excluded_type_counts.csv", index=False)

    source_counts = retained.dataset_id.value_counts(sort=False).sort_index()
    available_columns: Counter[str] = Counter()
    union: set[str] = set()
    for number, (dataset_id, n_stars) in enumerate(source_counts.items(), 1):
        symbols = source_symbols(dataset_id)
        union.update(symbols)
        available_columns.update(dict.fromkeys(symbols, int(n_stars)))
        if number % 100 == 0 or number == len(source_counts):
            print(f"GENE_AVAILABILITY {number}/{len(source_counts)} sources", flush=True)

    threshold = int(np.ceil(MIN_GENE_FRACTION * len(retained)))
    availability = pd.DataFrame(
        {
            "gene": sorted(union),
            "available_columns": [available_columns[g] for g in sorted(union)],
        }
    )
    availability["missing_columns"] = len(retained) - availability.available_columns
    availability["available_fraction"] = availability.available_columns / len(retained)
    availability["passes_90_percent"] = availability.available_columns.ge(threshold)
    availability.to_csv(OUTPUT / "gene_availability.csv", index=False)

    summary = {
        "complete": True,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "starting_supported_types": 71,
        "starting_supported_stars": 358_928,
        "starting_supported_maps": 16_366,
        "excluded_types": EXCLUDED_TYPES,
        "excluded_stars_before_map_filter": int(len(removed)),
        "retained_types_before_map_filter": int(retained.parent_id.nunique()),
        "retained_stars_before_map_filter": int(pre_map_filter_stars),
        "retained_maps_before_map_filter": int(pre_map_filter_maps),
        "minimum_stars_per_map": MIN_STARS_PER_MAP,
        "maps_dropped_after_type_exclusion": int(len(dropped_map_indices)),
        "retained_types": int(retained.parent_id.nunique()),
        "retained_stars": int(len(retained)),
        "retained_maps": int(retained.map_index.nunique()),
        "retained_sources": int(retained.dataset_id.nunique()),
        "retained_study_identifiers": int(retained.study_id.nunique()),
        "retained_donors": int(retained.donor_uid.nunique()),
        "union_unique_gene_symbols": int(len(union)),
        "minimum_gene_available_fraction": MIN_GENE_FRACTION,
        "minimum_gene_available_columns": threshold,
        "genes_passing_availability": int(availability.passes_90_percent.sum()),
        "universal_unique_genes": int(availability.available_columns.eq(len(retained)).sum()),
    }
    (OUTPUT / "cohort_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
