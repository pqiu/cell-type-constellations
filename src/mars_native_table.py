"""Strict reader for original MARS gene-by-cell integer-count tables."""
import numpy as np
from src.geo_inputs import read_gene_table


def read_native_table(raw, expected_genes, selected_cell_ids):
    lines = raw.splitlines()
    if len(lines) != len(expected_genes) + 1:
        raise ValueError('Native gene table row count changed')
    cells = lines[0].decode('utf-8-sig').split('\t')
    if not cells or any(not x for x in cells) or len(set(cells)) != len(cells):
        raise ValueError('Nonunique or empty original cell identifiers')
    for gene, line in zip(expected_genes, lines[1:]):
        name, sep, values = line.partition(b'\t')
        if not sep or name.decode() != gene or values.count(b'\t') != len(cells)-1:
            raise ValueError('Original native gene order or cell width changed')
        # Empty fields cannot silently collapse under whitespace parsing.
        if not values or values.startswith(b'\t') or values.endswith(b'\t') or b'\t\t' in values:
            raise ValueError('Missing source count')
    full, genes, barcodes = read_gene_table(raw)
    if genes[0].tolist() != expected_genes or barcodes.tolist() != cells:
        raise ValueError('Parsed native axis changed')
    if (full.data > np.iinfo(np.int32).max).any():
        raise ValueError('Source count exceeds exact int32 range')
    chosen = np.array([i for i, cell in enumerate(cells) if cell in selected_cell_ids], dtype=int)
    # All source values have passed the gate before selecting author-assigned cells.
    return full[chosen].astype(np.int32), [cells[i] for i in chosen], len(cells)
