"""Final source-processing entry point for c672834e-c3e3-49cb-81a5-4c844be4a975.

Downloads reviewed public inputs, applies the final SCimilarity annotation
contract, and writes sample-level cell-type stars.
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src.pbmc_pipeline import process
CONFIGS = [
    'configs/a199ca73-035d-44e2-9893-4c493151db21.json'
]

def main() -> None:
    for relative in CONFIGS:
        process(ROOT / relative)

if __name__ == "__main__":
    main()
