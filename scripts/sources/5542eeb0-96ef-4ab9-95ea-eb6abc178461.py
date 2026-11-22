"""Final source-processing entry point for 5542eeb0-96ef-4ab9-95ea-eb6abc178461.

Downloads reviewed public inputs, applies the final SCimilarity annotation
contract, and writes sample-level cell-type stars.
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src.pbmc_pipeline import process
CONFIGS = [
    'configs/4a5b00e0-1ba3-4fd4-af89-d3512eb20720.json',
    'configs/d9f3a82c-24bb-4277-b715-8ec92a33dd32.json'
]

def main() -> None:
    for relative in CONFIGS:
        process(ROOT / relative)

if __name__ == "__main__":
    main()
