"""Final source-processing entry point for E-MTAB-15719.

Downloads reviewed public inputs, applies the final SCimilarity annotation
contract, and writes sample-level cell-type stars.
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src.pbmc_pipeline import process
CONFIGS = [
    'configs/E-MTAB-15719_LC1-3.json',
    'configs/E-MTAB-15719_LC1_1-2.json',
    'configs/E-MTAB-15719_LC2.json',
    'configs/E-MTAB-15719_LC3.json'
]

def main() -> None:
    for relative in CONFIGS:
        process(ROOT / relative)

if __name__ == "__main__":
    main()
