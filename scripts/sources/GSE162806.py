"""Final source-processing entry point for GSE162806.

Downloads reviewed public inputs, applies the final SCimilarity annotation
contract, and writes sample-level cell-type stars.
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src.pbmc_pipeline import process
CONFIGS = [
    'configs/GSE162806_genes001.json',
    'configs/GSE162806_genes002.json',
    'configs/GSE162806_genes003.json',
    'configs/GSE162806_genes004.json',
    'configs/GSE162806_genes005.json',
    'configs/GSE162806_genes006.json',
    'configs/GSE162806_genes007.json',
    'configs/GSE162806_genes008.json',
    'configs/GSE162806_genes009.json',
    'configs/GSE162806_genes010.json'
]

def main() -> None:
    for relative in CONFIGS:
        process(ROOT / relative)

if __name__ == "__main__":
    main()
