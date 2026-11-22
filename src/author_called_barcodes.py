"""Select source-called droplets without substituting corrected expression."""
import numpy as np
import pandas as pd

def select_called_counts(raw_counts, raw_barcodes, called_barcodes):
    raw=pd.Index(raw_barcodes);called=pd.Index(called_barcodes)
    if not raw.is_unique or not called.is_unique or len(called)==0:
        raise ValueError('Original and author-called barcodes must be unique/nonempty')
    if len(raw)!=raw_counts.shape[0]:raise ValueError('Raw barcode/count dimensions differ')
    positions=raw.get_indexer(called)
    if (positions<0).any():raise ValueError('Author-called barcode absent from original counts')
    return raw_counts[positions].tocsr(), positions
