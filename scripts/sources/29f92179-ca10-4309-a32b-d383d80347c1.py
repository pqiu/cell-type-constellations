"""Final source-processing entry point for 29f92179-ca10-4309-a32b-d383d80347c1.

Downloads reviewed public inputs, applies the final SCimilarity annotation
contract, and writes sample-level cell-type stars.
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src.pbmc_pipeline import process
CONFIGS = [
    'configs/01c93cf6-b695-4e30-a26e-121ae8b16a9e.json',
    'configs/50eb1e23-b8d4-4f76-a184-44e5541fa05a.json',
    'configs/79ef1959-a6b4-4cac-82ca-30feaec48df1.json',
    'configs/84230ea4-998d-4aa8-8456-81dd54ce23af.json',
    'configs/881fe679-c6e0-45a3-9427-c4e81be6921f.json',
    'configs/8a554710-08bc-4005-87cd-da9675bdc2e7.json',
    'configs/c874f155-9bf9-4928-b821-f52c876b3e48.json',
    'configs/cab0bc48-744c-461b-aaf4-7bf2cb7af00d.json',
    'configs/db59611b-42de-4035-93aa-1ed39f38b467.json',
    'configs/ea786a06-5855-48b7-80d7-0313a21a2044.json',
    'configs/ed9e9f96-4f08-49d2-bef5-b2c29adf3edc.json'
]

def main() -> None:
    for relative in CONFIGS:
        process(ROOT / relative)

if __name__ == "__main__":
    main()
