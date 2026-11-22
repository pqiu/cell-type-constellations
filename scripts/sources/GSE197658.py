"""Final source-processing entry point for GSE197658.

Downloads reviewed public inputs, applies the final SCimilarity annotation
contract, and writes sample-level cell-type stars.
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src.pbmc_pipeline import process
CONFIGS = [
    'configs/GSE197658_RNA_18A.json',
    'configs/GSE197658_RNA_18D.json',
    'configs/GSE197658_RNA_18N.json',
    'configs/GSE197658_RNA_23A.json',
    'configs/GSE197658_RNA_23D.json',
    'configs/GSE197658_RNA_23N.json',
    'configs/GSE197658_RNA_30A.json',
    'configs/GSE197658_RNA_30D.json',
    'configs/GSE197658_RNA_30N.json',
    'configs/GSE197658_RNA_32A.json',
    'configs/GSE197658_RNA_32D.json',
    'configs/GSE197658_RNA_32N.json',
    'configs/GSE197658_RNA_35A.json',
    'configs/GSE197658_RNA_35D.json',
    'configs/GSE197658_RNA_35N.json',
    'configs/GSE197658_RNA_36A.json',
    'configs/GSE197658_RNA_36D.json',
    'configs/GSE197658_RNA_36N.json',
    'configs/GSE197658_RNA_37A.json',
    'configs/GSE197658_RNA_37D.json',
    'configs/GSE197658_RNA_37N.json',
    'configs/GSE197658_RNA_38A.json',
    'configs/GSE197658_RNA_38D.json',
    'configs/GSE197658_RNA_38N.json',
    'configs/GSE197658_RNA_39A.json',
    'configs/GSE197658_RNA_39D.json',
    'configs/GSE197658_RNA_39N.json',
    'configs/GSE197658_RNA_40A.json',
    'configs/GSE197658_RNA_40D.json',
    'configs/GSE197658_RNA_40N.json'
]

def main() -> None:
    for relative in CONFIGS:
        process(ROOT / relative)

if __name__ == "__main__":
    main()
