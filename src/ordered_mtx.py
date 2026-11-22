"""Bounded-memory reader for explicitly reviewed, row-ordered count matrices."""
import numpy as np
import pandas as pd


def ordered_count_chunks(stream, expected_shape, chunksize=1_000_000):
    header=stream.readline().strip()
    if header not in ('%%MatrixMarket matrix coordinate real general',
                      '%%MatrixMarket matrix coordinate integer general'):
        raise ValueError('Expected general coordinate count matrix')
    line=stream.readline()
    while line.startswith('%'):line=stream.readline()
    dims=tuple(map(int,line.split()))
    if len(dims)!=3 or dims[:2]!=tuple(expected_shape) or dims[2]<0:
        raise ValueError('Source dimensions differ from reviewed native axes')
    total=0;previous=-1;nrows,ncols=expected_shape
    for frame in pd.read_csv(stream,sep=r'\s+',header=None,names=['row','col','value'],
            dtype={'row':np.int64,'col':np.int64,'value':np.float64},chunksize=chunksize):
        row=frame.row.to_numpy()-1;col=frame.col.to_numpy()-1;value=frame.value.to_numpy()
        if (not np.isfinite(value).all() or (value<0).any()
                or not np.equal(value,np.floor(value)).all()
                or (value>np.iinfo(np.int32).max).any()):
            raise ValueError('Nonfinite, negative, fractional or out-of-range count')
        if (row<0).any() or (row>=nrows).any() or (col<0).any() or (col>=ncols).any():
            raise ValueError('Count coordinate outside native axes')
        keys=row*ncols+col
        if not len(keys) or keys[0]<=previous or (np.diff(keys)<=0).any():
            raise ValueError('Count coordinates must be strictly row-then-gene ordered, without duplicates')
        previous=int(keys[-1]);total+=len(row)
        yield row,col,value.astype(np.int32)
    if total!=dims[2]:raise ValueError('Declared nonzero-entry count differs from source records')
