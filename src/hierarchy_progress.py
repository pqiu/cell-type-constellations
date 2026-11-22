"""Durable per-source progress snapshot for parent-star regeneration."""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from scripts.inspect_pooled_star_axes import eligible_outputs


def write_progress(root: Path) -> dict:
    base = root / "results/hierarchy_star_maps"
    expected = {dataset_id: len(pd.read_csv(root / "results/star_maps" / dataset_id / "samples.csv",
                                           usecols=["sample_id"]))
                for dataset_id in eligible_outputs()}
    completed = []
    for path in base.glob("*/report.json"):
        if path.parent.name not in expected:
            continue
        report = json.loads(path.read_text(encoding="utf-8"))
        if int(report["samples"]) != expected[path.parent.name]:
            raise ValueError("Completed source has a different sample-view count than curation")
        completed.append({"dataset_id": path.parent.name,
                          "maps": int(report["samples"]),
                          "native_complete": "expression_native.h5" in report.get("output_sha256", {}),
                          "cells": int(report["direct_labeled_cells"])})
    failures = []
    for path in base.glob("*/failure.json"):
        if path.parent.name in expected and not any(row["dataset_id"] == path.parent.name and row["native_complete"]
                                                     for row in completed):
            failures.append(json.loads(path.read_text(encoding="utf-8")))
    native = [row for row in completed if row["native_complete"]]
    snapshot = {
        "updated_at_utc": datetime.now(timezone.utc).isoformat(),
        "expected_source_outputs": len(expected),
        "expected_maps": int(sum(expected.values())),
        "native_completed_outputs": len(native),
        "native_completed_maps": int(sum(row["maps"] for row in native)),
        "selected_gene_completed_outputs": len(completed),
        "selected_gene_completed_maps": int(sum(row["maps"] for row in completed)),
        "failed_outputs_current": len(failures),
        "completed": sorted(completed, key=lambda row: row["dataset_id"]),
        "failures": sorted(failures, key=lambda row: row["dataset_id"]),
    }
    base.mkdir(parents=True, exist_ok=True)
    (base / "progress.json").write_text(json.dumps(snapshot, indent=2) + "\n", encoding="utf-8")
    lines = ["# Hierarchy-pooled star-map progress", "",
             f"Updated {snapshot['updated_at_utc']}", "",
             f"- Full native-gene outputs: {len(native)}/{len(expected)}",
             f"- Full native-gene maps: {snapshot['native_completed_maps']}/{snapshot['expected_maps']}",
             f"- Source outputs requiring attention: {len(failures)}", "",
             "| Source output | Maps | Full native genes |", "| --- | ---: | :---: |"]
    lines.extend(f"| {row['dataset_id']} | {row['maps']} | {'yes' if row['native_complete'] else 'pending'} |"
                 for row in snapshot["completed"])
    if failures:
        lines += ["", "## Sources requiring attention", ""]
        lines.extend(f"- `{row['dataset_id']}`: {row['error_type']}: {row['error']}" for row in snapshot["failures"])
    (base / "progress.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return snapshot
