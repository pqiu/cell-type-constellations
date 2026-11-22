"""Final source-processing entry point for 77f9d7e9-5675-49c3-abed-ce02f39eef1b.

Downloads reviewed public inputs, applies the final SCimilarity annotation
contract, and writes sample-level cell-type stars.
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src.pbmc_pipeline import process
CONFIGS = [
    'configs/0a8e3443-c3e2-4918-84db-0495657d9175.json',
    'configs/0c0c49a2-b96b-432f-98a4-1decf9929149.json',
    'configs/3cb646ab-d411-4e5b-89eb-082a786c953b.json',
    'configs/420d0f5d-b7ba-4b1e-88d0-44e9f733febd.json',
    'configs/46104f0b-9af5-466a-ae0f-56b8dc1969a2.json',
    'configs/47e00f43-6176-44a5-b4d1-f88a1920fcab.json',
    'configs/5fde5c9c-5b1e-4df5-982a-3f8e7635161f.json',
    'configs/AIFI_human_immune_health_atlas_full.json',
    'configs/aa32098f-f7b0-43f2-a00a-601d190aae94.json'
]

def main() -> None:
    for relative in CONFIGS:
        process(ROOT / relative)

if __name__ == "__main__":
    main()
