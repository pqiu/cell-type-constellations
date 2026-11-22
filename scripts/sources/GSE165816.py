"""Final source-processing entry point for GSE165816.

Downloads reviewed public inputs, applies the final SCimilarity annotation
contract, and writes sample-level cell-type stars.
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src.pbmc_pipeline import process
CONFIGS = [
    'configs/GSE165816_G19.json',
    'configs/GSE165816_G20.json',
    'configs/GSE165816_G21.json',
    'configs/GSE165816_G22.json',
    'configs/GSE165816_G25.json',
    'configs/GSE165816_G26.json',
    'configs/GSE165816_G27.json',
    'configs/GSE165816_G30.json',
    'configs/GSE165816_G37.json',
    'configs/GSE165816_G47.json'
]

def main() -> None:
    for relative in CONFIGS:
        process(ROOT / relative)

if __name__ == "__main__":
    main()
