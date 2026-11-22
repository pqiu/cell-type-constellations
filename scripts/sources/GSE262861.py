"""Final source-processing entry point for GSE262861.

Downloads reviewed public inputs, applies the final SCimilarity annotation
contract, and writes sample-level cell-type stars.
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src.pbmc_pipeline import process
CONFIGS = [
    'configs/GSE262861_GSM8180643.json',
    'configs/GSE262861_GSM8180644.json',
    'configs/GSE262861_GSM8180645.json',
    'configs/GSE262861_GSM8180646.json',
    'configs/GSE262861_GSM8180647.json',
    'configs/GSE262861_GSM8180648.json',
    'configs/GSE262861_GSM8180649.json',
    'configs/GSE262861_GSM8180650.json',
    'configs/GSE262861_GSM8180651.json',
    'configs/GSE262861_GSM8180652.json',
    'configs/GSE262861_GSM8180653.json'
]

def main() -> None:
    for relative in CONFIGS:
        process(ROOT / relative)

if __name__ == "__main__":
    main()
