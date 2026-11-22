"""Final source-processing entry point for GSE239450.

Downloads reviewed public inputs, applies the final SCimilarity annotation
contract, and writes sample-level cell-type stars.
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src.pbmc_pipeline import process
CONFIGS = [
    'configs/GSM7665278_OB7-CTRL_GEX.json',
    'configs/GSM7665279_705385-SEV_GEX.json',
    'configs/GSM7665280_803763-ASX_GEX.json',
    'configs/GSM7665281_324058-ASX_GEX.json',
    'configs/GSM7665282_644394-CTRL_GEX.json',
    'configs/GSM7665283_915348-SEV_GEX.json',
    'configs/GSM7665284_729106-CTRL_GEX.json',
    'configs/GSM7665285_105199-ASX_GEX.json',
    'configs/GSM7665286_101607-SEV_GEX.json'
]

def main() -> None:
    for relative in CONFIGS:
        process(ROOT / relative)

if __name__ == "__main__":
    main()
