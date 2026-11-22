"""Shared release-level view selection for organ curation.

Dataset-level ``release_primary`` remains the default.  This manifest handles
the narrower case where selected views in an otherwise primary composite
source are superseded by better-provenanced original releases.
"""
from pathlib import Path

import pandas as pd


def release_view_exclusions(root: Path) -> set[tuple[str, str]]:
    path = root / "references/datasets/organ_curation_3000_v2/release_view_exclusions.csv"
    if not path.exists():
        return set()
    table = pd.read_csv(path, dtype=str).fillna("")
    required = {"dataset", "sample_id", "canonical_dataset", "canonical_sample_id", "reason", "evidence"}
    assert required.issubset(table.columns)
    keys = list(zip(table.dataset, table.sample_id))
    assert len(keys) == len(set(keys)), "Duplicate release-view exclusion"
    assert all(table.dataset.str.len()) and all(table.sample_id.str.len())
    return set(keys)


def exclude_release_views(table: pd.DataFrame, root: Path) -> pd.DataFrame:
    excluded = release_view_exclusions(root)
    if not excluded or table.empty:
        return table
    keep = [
        (str(dataset), str(sample_id)) not in excluded
        for dataset, sample_id in zip(table.dataset, table.sample_id)
    ]
    return table.loc[keep].copy()
