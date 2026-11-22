"""Final source-processing entry point for GSE285173.

Downloads reviewed public inputs, applies the final SCimilarity annotation
contract, and writes sample-level cell-type stars.
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src.pbmc_pipeline import process
CONFIGS = [
    'configs/GSE285173_GSM8697592.json',
    'configs/GSE285173_GSM8697593.json',
    'configs/GSE285173_GSM8697594.json',
    'configs/GSE285173_GSM8697595.json',
    'configs/GSE285173_GSM8697596.json',
    'configs/GSE285173_GSM8697597.json'
]

def main() -> None:
    for relative in CONFIGS:
        process(ROOT / relative)

if __name__ == "__main__":
    main()
