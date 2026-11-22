"""Build the broad 6-star Figure 3 on the finalized 66-star PBMC cohort."""
from pathlib import Path
import sys

import scripts.build_pbmc_6_star_broad_figure3 as figure


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "results/hierarchy_global_reanalysis/final_132/pbmc_66_within_map_gene_selection_v1"


def main() -> None:
    figure.BASE = BASE
    figure.RAW_DIR = BASE / "star_matrix_v1"
    figure.REF_DIR = BASE / "reference_map_v1"
    figure.OUT = BASE / "figure3_six_star_broad_v1"
    figure.RAW_MATRIX_NAME = "star_expression_2000_by_342559.float32.npy"
    figure.REFERENCE_MATRIX_NAME = "reference_expression_2000_by_66.float64.npy"
    if "--plot-only" in sys.argv:
        figure.plot_only()
    else:
        figure.main()


if __name__ == "__main__":
    main()
