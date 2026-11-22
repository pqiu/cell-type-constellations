"""Read wide original gene-by-cell tables without retaining decompressed text."""
import csv
import numpy as np
import pandas as pd
from scipy import sparse


def read_gene_table_stream(stream,separator='\t'):
    if separator not in ('\t',',',';',' '):raise ValueError('Unsupported reviewed gene-table separator')
    stream.seek(0)
    first=stream.readline();second=stream.readline();stream.seek(0)
    def decode(value):return value.decode('utf-8-sig') if isinstance(value,bytes) else value
    header=next(csv.reader([decode(first)],delimiter=separator,skipinitialspace=separator==' '))
    first_row=next(csv.reader([decode(second)],delimiter=separator,skipinitialspace=separator==' '))
    if len(first_row)==len(header)+1:barcodes=header
    elif len(first_row)==len(header):barcodes=header[1:]
    else:raise ValueError('Unexpected GEX table header width')
    if not barcodes or len(set(barcodes))!=len(barcodes) or any(not b for b in barcodes):
        raise ValueError('GEX table cell identifiers must be unique and nonempty')
    blocks=[];names=[]
    for frame in pd.read_csv(stream,sep=separator,index_col=0,chunksize=128,keep_default_na=False,
        skipinitialspace=separator==' ',encoding='utf-8-sig'):
        if frame.columns.tolist()!=barcodes:raise ValueError('GEX table cell order changed')
        values=frame.to_numpy(dtype=np.float64)
        if not np.isfinite(values).all() or (values<0).any() or not np.equal(values,np.floor(values)).all():
            raise ValueError('GEX table does not contain nonnegative integer counts')
        names.extend(frame.index.astype(str));blocks.append(sparse.csr_matrix(values))
    if not names or len(set(names))!=len(names) or any(not n.strip() for n in names):raise ValueError('Unique nonempty source gene identifiers required')
    return sparse.vstack(blocks,format='csr').T.tocsr(),pd.DataFrame({0:names,1:names}),pd.Series(barcodes)
