"""Build the 66-star pooled reference and one-shift-per-map aligned matrix."""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "results/hierarchy_global_reanalysis/final_132/pbmc_66_within_map_gene_selection_v1"
MATRIX_DIR = BASE / "star_matrix_v1"
REFERENCE_DIR = BASE / "reference_map_v1"
ALIGN_DIR = BASE / "whole_constellation_alignment_v1"
RAW_PATH = MATRIX_DIR / "star_expression_2000_by_342559.float32.npy"
ALIGNED_PATH = ALIGN_DIR / "aligned_expression_2000_by_342559.float32.npy"
SHIFTS_PATH = ALIGN_DIR / "map_shifts_2000_by_14230.float64.npy"


def sha256(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def atomic_json(path: Path, value: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def main() -> None:
    REFERENCE_DIR.mkdir(parents=True, exist_ok=True)
    ALIGN_DIR.mkdir(parents=True, exist_ok=True)
    raw = np.load(RAW_PATH, mmap_mode="r")
    stars = pd.read_csv(MATRIX_DIR / "star_columns.csv")
    genes = pd.read_csv(MATRIX_DIR / "gene_rows.csv").sort_values("matrix_row")
    if raw.shape != (2_000, 342_559) or stars.parent_id.nunique() != 66:
        raise ValueError("Unexpected finalized matrix inventory")
    if not np.array_equal(stars.matrix_column.to_numpy(), np.arange(len(stars))):
        raise ValueError("Star metadata are not in matrix-column order")

    types = (
        stars[["parent_id", "cell_type_name"]]
        .drop_duplicates()
        .sort_values(["parent_id", "cell_type_name"], kind="stable")
        .reset_index(drop=True)
    )
    types.insert(0, "reference_column", np.arange(len(types), dtype=np.int32))
    reference = np.empty((2_000, len(types)), dtype=np.float64)
    supports = []
    for row in types.itertuples(index=False):
        columns = stars.loc[stars.parent_id.eq(row.parent_id), "matrix_column"].to_numpy(dtype=np.int64)
        total = np.zeros(2_000, dtype=np.float64)
        for start in range(0, len(columns), 4_096):
            total += np.asarray(raw[:, columns[start:start + 4_096]], dtype=np.float64).sum(axis=1)
        reference[:, row.reference_column] = total / len(columns)
        supports.append(len(columns))
    types["observed_stars"] = supports
    types["observed_maps"] = supports
    types.to_csv(REFERENCE_DIR / "reference_star_columns.csv", index=False)
    np.save(REFERENCE_DIR / "reference_expression_2000_by_66.float64.npy", reference)
    table = pd.DataFrame(reference, columns=types.parent_id.tolist())
    table.insert(0, "gene", genes.gene.tolist())
    table.to_csv(REFERENCE_DIR / "reference_expression.csv.gz", index=False, compression="gzip")
    genes[["matrix_row", "gene"]].rename(columns={"matrix_row": "reference_row"}).to_csv(
        REFERENCE_DIR / "reference_gene_rows.csv", index=False
    )
    reference_receipt = {
        "complete": True,
        "verified": True,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "definition": "unweighted arithmetic mean across observed retained stars of each cell type",
        "shape": [2_000, 66],
        "genes": 2_000,
        "reference_stars": 66,
        "input_stars": len(stars),
        "input_maps": int(stars.map_index.nunique()),
        "minimum_type_support": int(min(supports)),
        "maximum_type_support": int(max(supports)),
        "all_values_finite": bool(np.isfinite(reference).all()),
        "input_matrix_sha256": sha256(RAW_PATH),
        "reference_matrix_sha256": sha256(
            REFERENCE_DIR / "reference_expression_2000_by_66.float64.npy"
        ),
    }
    atomic_json(REFERENCE_DIR / "completion.json", reference_receipt)

    maps = (
        stars.groupby("map_index", sort=True)
        .agg(
            sample_id=("sample_id", "first"), donor_uid=("donor_uid", "first"),
            dataset_id=("dataset_id", "first"), study_id=("study_id", "first"),
            first_matrix_column=("matrix_column", "min"),
            last_matrix_column=("matrix_column", "max"),
            observed_stars=("matrix_column", "size"),
        )
        .reset_index()
    )
    maps.insert(0, "alignment_map_column", np.arange(len(maps), dtype=np.int32))
    if len(maps) != 14_230 or maps.observed_stars.min() < 10:
        raise ValueError("Unexpected retained map inventory")
    if not np.array_equal(
        maps.last_matrix_column - maps.first_matrix_column + 1, maps.observed_stars
    ):
        raise ValueError("Map columns are not contiguous")
    lookup = dict(zip(types.parent_id, types.reference_column))
    aligned = np.lib.format.open_memmap(
        ALIGNED_PATH, mode="w+", dtype=np.float32, shape=raw.shape, fortran_order=True
    )
    shifts = np.lib.format.open_memmap(
        SHIFTS_PATH, mode="w+", dtype=np.float64, shape=(2_000, len(maps)), fortran_order=True
    )
    for completed, row in enumerate(maps.itertuples(index=False), 1):
        start = int(row.first_matrix_column)
        stop = int(row.last_matrix_column) + 1
        block = np.asarray(raw[:, start:stop], dtype=np.float64)
        ref_columns = np.fromiter(
            (lookup[x] for x in stars.loc[start:stop - 1, "parent_id"]),
            dtype=np.int64, count=stop - start,
        )
        shift = np.mean(block - reference[:, ref_columns], axis=1, dtype=np.float64)
        shifts[:, completed - 1] = shift
        aligned[:, start:stop] = block - shift[:, None]
        if completed % 1_000 == 0 or completed == len(maps):
            aligned.flush(); shifts.flush()
            print(f"ALIGNED {completed:,}/{len(maps):,} maps", flush=True)
    del aligned, shifts
    maps.to_csv(ALIGN_DIR / "map_columns.csv", index=False)
    (ALIGN_DIR / "star_columns.csv").write_bytes((MATRIX_DIR / "star_columns.csv").read_bytes())
    (ALIGN_DIR / "gene_rows.csv").write_bytes((MATRIX_DIR / "gene_rows.csv").read_bytes())

    aligned = np.load(ALIGNED_PATH, mmap_mode="r")
    maximum_error = 0.0
    for row in maps.itertuples(index=False):
        start = int(row.first_matrix_column)
        stop = int(row.last_matrix_column) + 1
        ref_columns = np.fromiter(
            (lookup[x] for x in stars.loc[start:stop - 1, "parent_id"]),
            dtype=np.int64, count=stop - start,
        )
        error = np.max(np.abs(np.mean(
            np.asarray(aligned[:, start:stop], dtype=np.float64) - reference[:, ref_columns], axis=1
        )))
        maximum_error = max(maximum_error, float(error))
    if not np.isfinite(aligned).all() or maximum_error > 2e-6:
        raise ValueError(f"Alignment verification failed: {maximum_error}")
    alignment_receipt = {
        "complete": True,
        "verified": True,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "definition": "one equal-star whole-constellation translation per retained map",
        "target_aware": False,
        "ancestor_descendant_exclusions": False,
        "iou_exclusions": False,
        "aligned_matrix_shape": list(aligned.shape),
        "shift_matrix_shape": [2_000, len(maps)],
        "maps": len(maps), "stars": len(stars), "genes": 2_000,
        "minimum_stars_per_map": int(maps.observed_stars.min()),
        "maximum_stars_per_map": int(maps.observed_stars.max()),
        "all_aligned_values_finite": True,
        "maximum_aligned_map_mean_displacement_error": maximum_error,
        "input_matrix_sha256": sha256(RAW_PATH),
        "reference_matrix_sha256": reference_receipt["reference_matrix_sha256"],
        "aligned_matrix_sha256": sha256(ALIGNED_PATH),
        "shift_matrix_sha256": sha256(SHIFTS_PATH),
    }
    atomic_json(ALIGN_DIR / "completion.json", alignment_receipt)
    print(json.dumps({"reference": reference_receipt, "alignment": alignment_receipt}, indent=2))


if __name__ == "__main__":
    main()
