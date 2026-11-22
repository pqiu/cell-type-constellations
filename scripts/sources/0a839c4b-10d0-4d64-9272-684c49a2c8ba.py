"""Final source-processing entry point for 0a839c4b-10d0-4d64-9272-684c49a2c8ba.

Downloads reviewed public inputs, applies the final SCimilarity annotation
contract, and writes sample-level cell-type stars.
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src.pbmc_pipeline import process
CONFIGS = [
    'configs/9dbab10c-118d-496b-966a-67f1763a6b7d.json',
    'configs/GSE158055_original_additional_views.json'
]

def main() -> None:
    for relative in CONFIGS:
        process(ROOT / relative)

if __name__ == "__main__":
    main()
