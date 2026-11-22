"""Final source-processing entry point for GSE266488.

Downloads reviewed public inputs, applies the final SCimilarity annotation
contract, and writes sample-level cell-type stars.
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src.pbmc_pipeline import process
CONFIGS = [
    'configs/GSE266488_GSM8178020.json',
    'configs/GSE266488_GSM8178021.json',
    'configs/GSE266488_GSM8178022.json',
    'configs/GSE266488_GSM8178023.json',
    'configs/GSE266488_GSM8178024.json',
    'configs/GSE266488_GSM8178025.json',
    'configs/GSE266488_GSM8178026.json',
    'configs/GSE266488_GSM8178027.json',
    'configs/GSE266488_GSM8178028.json',
    'configs/GSE266488_GSM8178029.json',
    'configs/GSE266488_GSM8178030.json',
    'configs/GSE266488_GSM8248050.json',
    'configs/GSE266488_GSM8248051.json',
    'configs/GSE266488_GSM8248052.json',
    'configs/GSE266488_GSM8248053.json',
    'configs/GSE266488_GSM8248054.json',
    'configs/GSE266488_GSM8248055.json',
    'configs/GSE266488_GSM8248056.json',
    'configs/GSE266488_GSM8248057.json',
    'configs/GSE266488_GSM8248058.json',
    'configs/GSE266488_GSM8248059.json',
    'configs/GSE266488_GSM8248060.json',
    'configs/GSE266488_GSM8248061.json',
    'configs/GSE266488_GSM8248062.json',
    'configs/GSE266488_GSM8248063.json',
    'configs/GSE266488_GSM8248064.json',
    'configs/GSE266488_GSM8248065.json',
    'configs/GSE266488_GSM8248066.json',
    'configs/GSE266488_GSM8248067.json',
    'configs/GSE266488_GSM8248068.json',
    'configs/GSE266488_GSM8248069.json',
    'configs/GSE266488_GSM8248070.json',
    'configs/GSE266488_GSM8248071.json',
    'configs/GSE266488_GSM8248072.json',
    'configs/GSE266488_GSM8248073.json',
    'configs/GSE266488_GSM8248074.json',
    'configs/GSE266488_GSM8248075.json',
    'configs/GSE266488_GSM8248076.json',
    'configs/GSE266488_GSM8248077.json',
    'configs/GSE266488_GSM8248078.json',
    'configs/GSE266488_GSM8248079.json',
    'configs/GSE266488_GSM8248080.json',
    'configs/GSE266488_GSM8248081.json',
    'configs/GSE266488_GSM8248082.json',
    'configs/GSE266488_GSM8248083.json',
    'configs/GSE266488_GSM8248084.json',
    'configs/GSE266488_GSM8248085.json',
    'configs/GSE266488_GSM8248086.json',
    'configs/GSE266488_GSM8248087.json',
    'configs/GSE266488_GSM8248088.json',
    'configs/GSE266488_GSM8248089.json',
    'configs/GSE266488_GSM8248090.json',
    'configs/GSE266488_GSM8248091.json',
    'configs/GSE266488_GSM8248092.json',
    'configs/GSE266488_GSM8248093.json',
    'configs/GSE266488_GSM8248094.json',
    'configs/GSE266488_GSM8248095.json',
    'configs/GSE266488_GSM8248096.json',
    'configs/GSE266488_GSM8248097.json'
]

def main() -> None:
    for relative in CONFIGS:
        process(ROOT / relative)

if __name__ == "__main__":
    main()
