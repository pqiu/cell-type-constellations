"""Final source-processing entry point for GSE316301.

Downloads reviewed public inputs, applies the final SCimilarity annotation
contract, and writes sample-level cell-type stars.
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src.pbmc_pipeline import process
CONFIGS = [
    'configs/GSE316301_GSM9449215.json',
    'configs/GSE316301_GSM9449216.json',
    'configs/GSE316301_GSM9449217.json',
    'configs/GSE316301_GSM9449218.json',
    'configs/GSE316301_GSM9449219.json',
    'configs/GSE316301_GSM9449220.json',
    'configs/GSE316301_GSM9449221.json',
    'configs/GSE316301_GSM9449222.json',
    'configs/GSE316301_GSM9449223.json',
    'configs/GSE316301_GSM9449224.json',
    'configs/GSE316301_GSM9449225.json',
    'configs/GSE316301_GSM9449226.json',
    'configs/GSE316301_GSM9449227.json',
    'configs/GSE316301_GSM9449228.json',
    'configs/GSE316301_GSM9449229.json',
    'configs/GSE316301_GSM9449230.json',
    'configs/GSE316301_GSM9449231.json',
    'configs/GSE316301_GSM9449232.json',
    'configs/GSE316301_GSM9449233.json',
    'configs/GSE316301_GSM9449234.json',
    'configs/GSE316301_GSM9449235.json',
    'configs/GSE316301_GSM9449236.json',
    'configs/GSE316301_GSM9449237.json',
    'configs/GSE316301_GSM9449238.json'
]

def main() -> None:
    for relative in CONFIGS:
        process(ROOT / relative)

if __name__ == "__main__":
    main()
