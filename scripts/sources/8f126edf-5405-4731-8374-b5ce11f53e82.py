"""Final source-processing entry point for 8f126edf-5405-4731-8374-b5ce11f53e82.

Downloads reviewed public inputs, applies the final SCimilarity annotation
contract, and writes sample-level cell-type stars.
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src.pbmc_pipeline import process
CONFIGS = [
    'configs/ebc2e1ff-c8f9-466a-acf4-9d291afaf8b3.json'
]

def main() -> None:
    for relative in CONFIGS:
        process(ROOT / relative)

if __name__ == "__main__":
    main()
