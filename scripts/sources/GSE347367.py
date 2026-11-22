"""Final source-processing entry point for GSE347367.

Downloads reviewed public inputs, applies the final SCimilarity annotation
contract, and writes sample-level cell-type stars.
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src.pbmc_pipeline import process
CONFIGS = [
    'configs/GSE347367_BNT_V07_ID0.json',
    'configs/GSE347367_BNT_V07_ID1.json',
    'configs/GSE347367_BNT_V09_ID0.json',
    'configs/GSE347367_BNT_V09_ID1.json',
    'configs/GSE347367_BNT_V30_ID0.json',
    'configs/GSE347367_BNT_V30_ID1.json',
    'configs/GSE347367_BNT_V30_ID2.json',
    'configs/GSE347367_BNT_V43_ID0.json',
    'configs/GSE347367_BNT_V43_ID1.json',
    'configs/GSE347367_BNT_V43_ID3.json',
    'configs/GSE347367_BNT_V46_ID0.json',
    'configs/GSE347367_BNT_V46_ID1.json',
    'configs/GSE347367_BNT_V46_ID3.json',
    'configs/GSE347367_BNT_V69_ID0.json',
    'configs/GSE347367_BNT_V69_ID1.json',
    'configs/GSE347367_BNT_V69_ID2.json',
    'configs/GSE347367_BNT_V69_ID3.json',
    'configs/GSE347367_BNT_V74_ID0.json',
    'configs/GSE347367_BNT_V74_ID1.json',
    'configs/GSE347367_BNT_V74_ID3.json',
    'configs/GSE347367_Moderna_V01_ID1.json',
    'configs/GSE347367_Moderna_V02_ID0.json',
    'configs/GSE347367_Moderna_V10_ID1.json',
    'configs/GSE347367_Moderna_V10_ID2.json',
    'configs/GSE347367_Moderna_V10_ID3.json',
    'configs/GSE347367_Moderna_V14_ID1.json',
    'configs/GSE347367_Moderna_V16_ID0.json',
    'configs/GSE347367_Moderna_V16_ID1.json',
    'configs/GSE347367_Moderna_V20_ID0.json',
    'configs/GSE347367_Moderna_V20_ID1.json',
    'configs/GSE347367_Moderna_V23_ID0.json',
    'configs/GSE347367_Moderna_V23_ID1.json',
    'configs/GSE347367_Moderna_V23_ID3.json',
    'configs/GSE347367_Moderna_V24_ID0.json',
    'configs/GSE347367_Moderna_V24_ID1.json',
    'configs/GSE347367_Moderna_V26_ID0.json',
    'configs/GSE347367_Moderna_V26_ID1.json',
    'configs/GSE347367_Moderna_V26_ID3.json',
    'configs/GSE347367_Moderna_V27_ID1.json',
    'configs/GSE347367_Moderna_V36_ID0.json',
    'configs/GSE347367_Moderna_V36_ID1.json',
    'configs/GSE347367_Moderna_V38_ID1.json',
    'configs/GSE347367_Moderna_V41_ID1.json',
    'configs/GSE347367_Moderna_V44_ID0.json',
    'configs/GSE347367_Moderna_V44_ID1.json',
    'configs/GSE347367_Moderna_V44_ID3.json',
    'configs/GSE347367_Moderna_V45_ID0.json',
    'configs/GSE347367_Moderna_V45_ID1.json',
    'configs/GSE347367_Moderna_V45_ID3.json',
    'configs/GSE347367_Moderna_V47_ID0.json',
    'configs/GSE347367_Moderna_V47_ID1.json',
    'configs/GSE347367_Moderna_V47_ID2.json',
    'configs/GSE347367_Moderna_V50_ID0.json',
    'configs/GSE347367_Moderna_V50_ID1.json',
    'configs/GSE347367_Moderna_V50_ID2.json',
    'configs/GSE347367_Moderna_V51_ID0.json',
    'configs/GSE347367_Moderna_V52_ID0.json',
    'configs/GSE347367_Moderna_V52_ID1.json',
    'configs/GSE347367_Moderna_V53_ID0.json',
    'configs/GSE347367_Moderna_V53_ID1.json',
    'configs/GSE347367_Moderna_V53_ID3.json',
    'configs/GSE347367_Moderna_V54_ID0.json',
    'configs/GSE347367_Moderna_V54_ID1.json',
    'configs/GSE347367_Moderna_V54_ID2.json',
    'configs/GSE347367_Moderna_V54_ID3.json',
    'configs/GSE347367_Moderna_V55_ID1.json',
    'configs/GSE347367_Moderna_V56_ID0.json',
    'configs/GSE347367_Moderna_V56_ID1.json',
    'configs/GSE347367_Moderna_V57_ID0.json',
    'configs/GSE347367_Moderna_V57_ID1.json',
    'configs/GSE347367_Moderna_V60_ID1.json',
    'configs/GSE347367_Moderna_V60_ID2.json',
    'configs/GSE347367_Moderna_V60_ID3.json',
    'configs/GSE347367_Moderna_V61_ID1.json',
    'configs/GSE347367_Moderna_V61_ID3.json',
    'configs/GSE347367_Moderna_V65_ID1.json',
    'configs/GSE347367_Moderna_V77_ID1.json'
]

def main() -> None:
    for relative in CONFIGS:
        process(ROOT / relative)

if __name__ == "__main__":
    main()
