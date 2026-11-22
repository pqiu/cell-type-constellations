"""Exact observed-cell joins for separate integer sample-tag CSV tables."""
import re
import numpy as np
import pandas as pd


def read_observed_tag_csv(path,rna_barcodes):
    table=pd.read_csv(path,dtype=str,keep_default_na=False,encoding='utf-8-sig')
    if not {'cell_barcode','total'}.issubset(table):raise ValueError('Expected explicit cell_barcode and total fields')
    barcode=table.pop('cell_barcode')
    if not barcode.is_unique or not barcode.str.fullmatch('[ACGT]+').all():raise ValueError('Unique original bare DNA barcodes required')
    table.index=pd.Index(barcode,name='source_tag_cell_barcode')
    columns=[c for c in table if c!='total']
    if len(columns)<2 or not all(re.fullmatch(r'(?:PHI|HV)-\d+_.+',c) for c in columns):
        raise ValueError('Explicit original clinical sample-tag labels required')
    for c in table:
        if not table[c].str.fullmatch(r'\d+').all():raise ValueError('Native nonnegative integer tag counts required')
    values=table.astype(np.int64)
    if (values<0).any().any() or not np.array_equal(values[columns].sum(axis=1),values['total']):
        raise ValueError('Reported total must match all native tag columns')
    rna=pd.Index(rna_barcodes.astype(str))
    if not rna.is_unique or not rna.str.fullmatch(r'[ACGT]+-1').all():raise ValueError('Unique RNA barcodes with explicit -1 suffix required')
    bare=rna.str.slice(stop=-2)
    observed=bare.isin(values.index)
    selected=values.loc[bare[observed],columns].copy();selected.index=rna[observed]
    review=dict(source_tag_cells=len(values),source_RNA_cells=len(rna),RNA_cells_with_observed_tag_rows=int(observed.sum()),
        RNA_cells_without_observed_tag_rows=int((~observed).sum()),tag_rows_without_RNA_barcodes=int((~values.index.isin(bare)).sum()),
        barcode_join='Exact original DNA barcode after explicit removal of verified RNA -1 suffix; no row-order mapping or missing-row zero fill.',
        clinical_alias_policy='Clinical PHI/HV tag labels remain original opaque source codes. No mapping to Subject identifiers is inferred from tag or column order.')
    return selected,observed,review
