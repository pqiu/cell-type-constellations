"""Read uncorrected Cell Ranger RNA counts for exact source-called barcodes."""
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

def read_selected(path,barcodes,library_id):
    selected=pd.Index(barcodes).astype(str)
    if not selected.is_unique:
        raise ValueError('Repeated called-cell barcode')
    with h5py.File(path,'r') as f:
        if set(f)!= {'matrix'} or not str(f.attrs.get('software_version','')).startswith('cellranger-'):
            raise ValueError('Expected original Cell Ranger matrix, not corrected count output')
        declared=[value.decode() if isinstance(value,bytes) else str(value) for value in f.attrs['library_ids']]
        if declared!=[library_id]:
            raise ValueError('Cell Ranger source library differs')
        group=f['matrix']; features=group['features']
        source_barcodes=pd.Index(group['barcodes'].asstr()[:])
        if not source_barcodes.is_unique:
            raise ValueError('Repeated raw matrix barcode')
        columns=source_barcodes.get_indexer(selected)
        if (columns<0).any():
            raise ValueError('Called cell missing from original matrix')
        raw=sparse.csc_matrix((group['data'][:],group['indices'][:],group['indptr'][:]),shape=tuple(group['shape'][:]))
        if not raw.has_canonical_format:
            raise ValueError('Invalid original sparse coordinates')
        rna=features['feature_type'].asstr()[:]=='Gene Expression'
        genes=pd.Index(features['id'].asstr()[:][rna])
        symbols=features['name'].asstr()[:][rna]
        if not genes.is_unique or not len(genes):
            raise ValueError('Invalid original RNA gene identities')
        counts=raw[rna][:,columns].T.tocsr()
        counts.eliminate_zeros()
        if not np.isfinite(counts.data).all() or (counts.data<0).any() or not np.equal(counts.data,np.floor(counts.data)).all():
            raise ValueError('Original RNA values are not counts')
        return counts,genes,symbols,dict(raw_barcodes=len(source_barcodes),called_barcodes=len(selected),source_software=str(f.attrs['software_version']))
