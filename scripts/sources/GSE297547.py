"""Final source-processing entry point for GSE297547.

Downloads reviewed public inputs, applies the final SCimilarity annotation
contract, and writes sample-level cell-type stars.
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src.pbmc_pipeline import process
CONFIGS = [
    'configs/GSE297547_C001.json',
    'configs/GSE297547_C002.json',
    'configs/GSE297547_C003.json',
    'configs/GSE297547_C004.json',
    'configs/GSE297547_C005.json',
    'configs/GSE297547_C006.json',
    'configs/GSE297547_C007.json',
    'configs/GSE297547_C008.json',
    'configs/GSE297547_C009.json',
    'configs/GSE297547_C010.json',
    'configs/GSE297547_MG001.json',
    'configs/GSE297547_MG002.json',
    'configs/GSE297547_MG003.json',
    'configs/GSE297547_MG004.json',
    'configs/GSE297547_MG005.json',
    'configs/GSE297547_MG006.json',
    'configs/GSE297547_MG007.json',
    'configs/GSE297547_MG008.json',
    'configs/GSE297547_MG009.json'
]

def main() -> None:
    for relative in CONFIGS:
        process(ROOT / relative)

if __name__ == "__main__":
    main()
