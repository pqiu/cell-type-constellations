"""Final source-processing entry point for GSE295527.

Downloads reviewed public inputs, applies the final SCimilarity annotation
contract, and writes sample-level cell-type stars.
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src.pbmc_pipeline import process
CONFIGS = [
    'configs/GSE295527_archive_genes001.json',
    'configs/GSE295527_archive_genes002.json',
    'configs/GSE295527_flex_gsm8951957.json',
    'configs/GSE295527_flex_gsm8951958.json',
    'configs/GSE295527_flex_gsm8951959.json',
    'configs/GSE295527_flex_gsm8951960.json',
    'configs/GSE295527_genes001.json'
]

def main() -> None:
    for relative in CONFIGS:
        process(ROOT / relative)

if __name__ == "__main__":
    main()
