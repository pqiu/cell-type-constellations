"""Final source-processing entry point for GSE132065.

Downloads reviewed public inputs, applies the final SCimilarity annotation
contract, and writes sample-level cell-type stars.
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src.pbmc_pipeline import process
CONFIGS = [
    'configs/GSE132065_clinical_RT.json',
    'configs/GSE132065_healthy_JULIA03.json'
]

def main() -> None:
    for relative in CONFIGS:
        process(ROOT / relative)

if __name__ == "__main__":
    main()
