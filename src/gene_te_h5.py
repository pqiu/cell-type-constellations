"""Select original Ensembl gene rows from an explicitly reviewed gene/TE reference.

RepeatMasker rows are RNA measurements, not additional genes. The selected
Ensembl gene universe defines the gene-expression library and centroid axes.
Both gene and excluded TE library totals remain recorded for every barcode.
"""
import hashlib
import json
import re

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

TE_PREFIX = 'hg19_TE_repmask_LTRm_s_20140131#'


def read_matrix(matrix):
    f=matrix['features']; n=len(f['id']); gene_ids=[]; te_count=0; entered_te=False
    for start in range(0,n,100000):
        stop=min(n,start+100000)
        ids=f['id'].asstr()[start:stop]
        if not np.all(f['feature_type'].asstr()[start:stop]=='Gene Expression'):
            raise ValueError('Unexpected feature type in reviewed gene/TE RNA reference')
        if not np.all(f['genome'].asstr()[start:stop]=='hg19_genome'):
            raise ValueError('Unexpected gene/TE source genome')
        for value in ids:
            if re.fullmatch(r'ENSG[0-9]+(?:\.[0-9]+)?',value):
                if entered_te:raise ValueError('Gene rows must precede the reviewed TE rows')
                gene_ids.append(value)
            else:
                entered_te=True
                if value!=TE_PREFIX+str(te_count):raise ValueError('Unrecognized or reordered non-gene feature')
                te_count+=1
    if not gene_ids or not te_count or len(set(gene_ids))!=len(gene_ids):
        raise ValueError('Unique genes and explicit sequential RepeatMasker rows required')
    names=f['name'].asstr()[:len(gene_ids)]
    if any(not str(v).strip() for v in names):raise ValueError('Empty native gene symbol')
    bars=pd.Series(matrix['barcodes'].asstr()[()])
    if not bars.is_unique or bars.str.strip().eq('').any():raise ValueError('Ambiguous original barcodes')
    values=matrix['data'][()]
    if not np.isfinite(values).all() or (values<0).any() or not np.equal(values,np.floor(values)).all():
        raise ValueError('Invalid original gene or TE RNA counts')
    all_counts=sparse.csc_matrix((values,matrix['indices'][()],matrix['indptr'][()]),shape=tuple(matrix['shape'][()])).T.tocsr()
    all_counts.check_format(full_check=True)
    if all_counts.shape!=(len(bars),n):raise ValueError('Original feature/barcode dimensions disagree')
    all_counts.sum_duplicates();all_counts.sort_indices();all_counts.eliminate_zeros()
    counts=all_counts[:,:len(gene_ids)].tocsr()
    total=np.asarray(all_counts.sum(axis=1,dtype=np.int64)).ravel()
    gene_total=np.asarray(counts.sum(axis=1,dtype=np.int64)).ravel()
    coordinate=hashlib.sha256(json.dumps([gene_ids,names.tolist()],separators=(',',':')).encode()).hexdigest()
    review=dict(rule='ensembl_genes_excluding_explicit_repeatmasker_TE',source_genome='hg19_genome',
        source_features=n,selected_Ensembl_gene_features=len(gene_ids),excluded_TE_features=te_count,
        original_all_feature_RNA_counts=int(total.sum()),selected_gene_RNA_counts=int(gene_total.sum()),
        excluded_TE_RNA_counts=int((total-gene_total).sum()),gene_coordinate_sha256=coordinate,
        normalization_denominator='Selected original Ensembl gene counts per cell; excluded repeat-element RNA does not enter gene library, MT fraction, annotation or star means.')
    libraries=pd.DataFrame(dict(source_original_all_feature_RNA_counts=total,source_excluded_TE_RNA_counts=total-gene_total))
    return counts,pd.DataFrame({0:gene_ids,1:names}),bars,review,libraries


def read(path):
    with h5py.File(path) as h:
        if list(h)!=['matrix']:raise ValueError('Unexpected source H5 layout')
        return read_matrix(h['matrix'])
