"""Final source-processing entry point for GSE317309.

Downloads reviewed public inputs, applies the final SCimilarity annotation
contract, and writes sample-level cell-type stars.
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src.pbmc_pipeline import process
CONFIGS = [
    'configs/GSE317309_d0_genes001.json',
    'configs/GSE317309_d21_genes001.json',
    'configs/GSE317309_d42_genes001.json',
    'configs/GSE317309_d63_genes001.json'
]

def main() -> None:
    for relative in CONFIGS:
        process(ROOT / relative)

if __name__ == "__main__":
    main()
