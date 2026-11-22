"""Final source-processing entry point for GSE190992.

Downloads reviewed public inputs, applies the final SCimilarity annotation
contract, and writes sample-level cell-type stars.
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src.pbmc_pipeline import process
CONFIGS = [
    'configs/GSE190992_native_ptid2.json',
    'configs/GSE190992_native_ptid4.json',
    'configs/GSE190992_native_ptid5.json',
    'configs/GSE190992_native_ptid6.json'
]

def main() -> None:
    for relative in CONFIGS:
        process(ROOT / relative)

if __name__ == "__main__":
    main()
