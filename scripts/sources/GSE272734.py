"""Final source-processing entry point for GSE272734.

Downloads reviewed public inputs, applies the final SCimilarity annotation
contract, and writes sample-level cell-type stars.
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src.pbmc_pipeline import process
CONFIGS = [
    'configs/GSE272734_genes001.json',
    'configs/GSE272734_genes002.json',
    'configs/GSE272734_genes003.json',
    'configs/GSE272734_genes004.json',
    'configs/GSE272734_genes005.json',
    'configs/GSE272734_genes006.json',
    'configs/GSE272734_genes007.json',
    'configs/GSE272734_genes008.json',
    'configs/GSE272734_genes009.json',
    'configs/GSE272734_genes010.json',
    'configs/GSE272734_genes011.json',
    'configs/GSE272734_genes012.json',
    'configs/GSE272734_genes013.json'
]

def main() -> None:
    for relative in CONFIGS:
        process(ROOT / relative)

if __name__ == "__main__":
    main()
