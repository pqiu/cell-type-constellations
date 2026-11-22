"""Read native gene-by-cell tables with repeated gene metadata between cell blocks."""
import csv
import io
import re

import numpy as np
import pandas as pd
from scipy import sparse


def read_table(raw):
    """Preserve literal gene IDs and symbols, verifying every repeated metadata field."""
    text = raw.decode('utf-8')
    reader = csv.reader(io.StringIO(text), delimiter='\t')
    header = next(reader)
    if header[:3] != ['Gene ID_1', 'Gene ID_2', 'Gene Symbol']:
        raise ValueError('Unexpected reviewed gene metadata columns')
    metadata_positions = [i for i, name in enumerate(header) if name in {'Gene ID_1', 'Gene ID_2', 'Gene Symbol'}]
    if any(header[i] == 'Gene ID_1' for i in metadata_positions if i != 0):
        raise ValueError('Unexpected repeated primary gene ID column')
    cell_positions = [i for i in range(len(header)) if i not in metadata_positions]
    barcodes = [header[i] for i in cell_positions]
    if not barcodes or len(barcodes) != len(set(barcodes)) or any(not b.strip() for b in barcodes):
        raise ValueError('Cell identifiers must be unique and nonempty')
    gene_rows = []; blocks = []
    boundaries = metadata_positions + [len(header)]
    segments = [(a + 1, b) for a, b in zip(boundaries, boundaries[1:]) if b > a + 1]
    buffer = np.empty((128, len(barcodes)), dtype=np.float64); used = 0
    for row in reader:
        if len(row) != len(header):
            raise ValueError('Expression row width differs from source header')
        for i in metadata_positions:
            original = {'Gene ID_1': 0, 'Gene ID_2': 1, 'Gene Symbol': 2}[header[i]]
            if row[i] != row[original]:
                raise ValueError('Repeated gene metadata differs from the original column')
        fields = []
        for a, b in segments: fields.extend(row[a:b])
        numeric_text = '\t'.join(fields)
        # These deposited tables encode counts as decimal integer tokens. Check
        # the entire string before fromstring, which otherwise accepts prefixes.
        if re.fullmatch(r'[0-9]+(?:\t[0-9]+)*', numeric_text) is None:
            raise ValueError('Native expression values must be finite nonnegative integer tokens')
        values = np.fromstring(numeric_text, sep='\t')
        if len(values) != len(barcodes) or not np.isfinite(values).all():
            raise ValueError('Invalid native expression row')
        if (values >= 2**53).any():
            raise ValueError('Count exceeds exact float64 integer range')
        gene_rows.append((row[0], row[2], row[1])); buffer[used] = values; used += 1
        if used == len(buffer):
            blocks.append(sparse.csr_matrix(buffer)); used = 0
    if used: blocks.append(sparse.csr_matrix(buffer[:used]))
    if not gene_rows:
        raise ValueError('Empty expression table')
    genes = pd.DataFrame(gene_rows, columns=[0, 1, 'source_gene_secondary_id'])
    if not genes[0].is_unique or genes[0].str.strip().eq('').any():
        raise ValueError('Native primary gene IDs must be unique and nonempty')
    # Missing source symbols do not invalidate a native Ensembl coordinate.
    # Keep its counts in native libraries/geometry; annotation ignores the empty symbol.
    counts = sparse.vstack(blocks).T.tocsr()
    return counts, genes, pd.Series(barcodes)
