"""Inspect common, unambiguous gene symbols in validated star-map exports."""
import collections
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STATUS = ROOT / "references/datasets/processing_status.csv"
CHECKPOINT = ROOT / "references/datasets/current_independent_validation_checkpoint.json"


def eligible_outputs():
    with STATUS.open(newline="", encoding="utf-8") as handle:
        ids = {r["dataset_id"] for r in csv.DictReader(handle) if r["status"] == "validated"}
    checkpoint = json.loads(CHECKPOINT.read_text(encoding="utf-8"))
    for item in checkpoint["independently_validated_additions"]:
        receipt = json.loads((ROOT / item["path"]).read_text(encoding="utf-8"))
        ids.update(Path(o["path"]).name for o in receipt["validated_outputs"])
    if len(ids) != checkpoint["current_independent_counts"]["processed_dataset_outputs"]:
        raise ValueError("Eligible output count does not match independent checkpoint")
    return sorted(ids)


def symbol_column(fields):
    for field in ("gene_symbols", "gene_symbol", "feature_name"):
        if field in fields:
            return field
    return "source_gene_id" if "source_gene_id" in fields else None


def main():
    ids = eligible_outputs()
    available = collections.Counter()
    unique = collections.Counter()
    no_symbols = []
    sizes = []
    for dataset_id in ids:
        path = ROOT / "results/star_maps" / dataset_id / "genes.csv"
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            field = symbol_column(reader.fieldnames)
            if field is None:
                no_symbols.append(dataset_id)
                continue
            counts = collections.Counter(r[field].strip() for r in reader if r[field] and r[field].strip())
        sizes.append(len(counts))
        available.update(counts.keys())
        unique.update(k for k, n in counts.items() if n == 1)
    n = len(ids) - len(no_symbols)
    print("eligible outputs", len(ids), "with symbols", n, "without", len(no_symbols))
    print("no symbol IDs", no_symbols)
    print("symbol counts min/median/max", min(sizes), sorted(sizes)[len(sizes)//2], max(sizes))
    for threshold in (1, .99, .95, .9, .8):
        print("fraction", threshold, "present", sum(v >= threshold*n for v in available.values()),
              "unique", sum(v >= threshold*n for v in unique.values()))


if __name__ == "__main__":
    main()
