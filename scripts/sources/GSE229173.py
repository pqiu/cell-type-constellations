"""Final source-processing entry point for GSE229173.

Downloads reviewed public inputs, applies the final SCimilarity annotation
contract, and writes sample-level cell-type stars.
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src.pbmc_pipeline import process
CONFIGS = [
    'configs/GSE229173_genes001.json',
    'configs/GSE229173_genes002.json',
    'configs/GSE229173_genes003.json',
    'configs/GSE229173_genes004.json',
    'configs/GSE229173_genes005.json',
    'configs/GSE229173_genes006.json'
]

def main() -> None:
    for relative in CONFIGS:
        process(ROOT / relative)

if __name__ == "__main__":
    main()
