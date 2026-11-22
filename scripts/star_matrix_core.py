"""Gene-axis helpers shared by final matrix construction."""
import hashlib
from pathlib import Path
import numpy as np
import pandas as pd
from scripts.inspect_pooled_star_axes import symbol_column

def digest(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()

def symbol_groups(path: Path) -> dict[str, np.ndarray]:
    genes = pd.read_csv(path, keep_default_na=False)
    field = symbol_column(genes.columns)
    names = genes[field].astype(str).str.strip().to_numpy()
    grouped = {}
    for row, name in enumerate(names):
        if name:
            grouped.setdefault(name, []).append(row)
    return {name: np.asarray(rows, dtype=np.int64) for name, rows in grouped.items()}
