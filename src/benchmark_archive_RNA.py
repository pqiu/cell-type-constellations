"""Read unchanged benchmark RNA archives; append native CSR libraries to H5AD."""
import gzip
import hashlib
import io
import json
import shutil
import tarfile
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
from threadpoolctl import threadpool_limits

from scripts.download_pbmc_asset import ROOT, LIMIT, footprint, download


def digest(path):
    with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


class CountSemanticsError(ValueError):
    def __init__(self,values,shape):
        finite=np.isfinite(values);fractional=finite & np.not_equal(values,np.floor(values))
        self.review=dict(source_matrix_shape=list(shape),entries=len(values),nonfinite_entries=int((~finite).sum()),
            negative_entries=int((values<0).sum()),fractional_entries=int(fractional.sum()),fractional_examples=values[fractional][:12].tolist(),
            interpretation='Native vendor matrix is not an integer UMI matrix. Not inferred normalized; no rounding/inverse transformation. Count semantics must be reviewed before annotation.')
        super().__init__('Original count semantics require review: '+json.dumps(self.review))


def read_capture(path, kind):
    if kind not in {'Parse','Fluent','Scale'}:raise ValueError('Explicit reviewed archive platform required')
    found={};members={}
    with tarfile.open(path,'r|*') as tar:
        for member in tar:
            if not member.isfile():continue
            name=member.name.rsplit('/',1)[-1]
            role=({'count_matrix.mtx':'matrix','all_genes.csv':'genes','cell_metadata.csv':'cells'} if kind=='Parse' else
                  {'matrix.mtx.gz':'matrix','features.tsv.gz':'genes','barcodes.tsv.gz':'cells'}).get(name)
            if role is None:continue
            if role in found:raise ValueError('Duplicate original count/coordinate member')
            members[role]=dict(name=member.name,bytes=member.size)
            with tar.extractfile(member) as raw:
                stream=gzip.GzipFile(fileobj=raw) if name.endswith('.gz') else raw
                if role=='matrix':
                    with threadpool_limits(8):found[role]=mmread(stream,spmatrix=True).tocsr()
                else:
                    data=stream.read();members[role]['decoded_sha256']=hashlib.sha256(data).hexdigest()
                    found[role]=pd.read_csv(io.BytesIO(data),sep=',' if kind=='Parse' else '\t',header=0 if kind=='Parse' else None,
                                            dtype=str,keep_default_na=False)
    if set(found)!={'matrix','genes','cells'}:raise ValueError('Missing original count/coordinate members')
    native=found['genes'];cells=found['cells'];x=found['matrix']
    x.check_format(full_check=True)
    if not np.isfinite(x.data).all() or (x.data<0).any() or not np.equal(x.data,np.floor(x.data)).all():raise CountSemanticsError(x.data,x.shape)
    if x.data.size and x.data.max()>=2**63:raise ValueError('Counts exceed exact signed integer storage')
    x.sum_duplicates();x.sort_indices();x.eliminate_zeros()
    if kind=='Parse':
        if native.columns.tolist()!=['gene_id','gene_name','genome'] or set(native.genome)!={'hg38'}:raise ValueError('Unreviewed Parse gene coordinates')
        if not {'bc_wells','sample','species','gene_count','tscp_count'}.issubset(cells) or set(cells.species)!={'hg38'}:raise ValueError('Invalid Parse cell metadata')
        if x.shape!=(len(cells),len(native)):raise ValueError('Original cells-by-genes orientation differs')
        np.testing.assert_array_equal(np.asarray(x.sum(axis=1)).ravel(),pd.to_numeric(cells.tscp_count))
        np.testing.assert_array_equal(np.diff(x.indptr),pd.to_numeric(cells.gene_count))
        ids=native.gene_id;symbols=native.gene_name;bars=cells.bc_wells
    else:
        if native.shape[1]!=3 or set(native[2])!={'Gene Expression'} or cells.shape[1]!=1:raise ValueError('Unreviewed triplet feature or barcode schema')
        if x.shape!=(len(native),len(cells)):raise ValueError('Original gene-by-cell dimensions differ')
        x=x.T.tocsr();x.sort_indices();ids=native[0];symbols=native[1];bars=cells[0]
        cells=pd.DataFrame({'barcode':bars})
    if not ids.is_unique or not bars.is_unique or ids.str.strip().eq('').any() or bars.str.strip().eq('').any():raise ValueError('Ambiguous original cell/gene coordinates')
    var=pd.DataFrame(dict(gene_symbols=symbols.to_numpy()),index=ids.to_numpy())
    axis_sha=hashlib.sha256(pd.DataFrame(dict(id=ids.to_numpy(),symbol=symbols.to_numpy())).to_csv(index=False).encode()).hexdigest()
    totals=np.asarray(x.sum(axis=1)).ravel();detected=np.diff(x.indptr)
    mito=np.asarray(x[:,symbols.str.startswith('MT-').to_numpy()].sum(axis=1)).ravel()
    qc=pd.DataFrame(dict(source_barcode=bars.to_numpy(),RNA_total=totals,RNA_detected_genes=detected,mito_RNA_total=mito))
    qc['RNA_QC_retained']=(totals>=500)&(detected>=200)&(mito/np.maximum(totals,1)<=.2)
    review=dict(platform=kind,source_cells=len(bars),source_genes=len(ids),native_mito_genes=int(symbols.str.startswith('MT-').sum()),
                nonzeros=int(x.nnz),total_RNA_counts=int(totals.sum()),gene_coordinate_sha256=axis_sha,members=members,
                integer_nonnegative_counts=True,source_native_symbol_blanks=int(symbols.eq('').sum()),RNA_QC_cells=int(qc.RNA_QC_retained.sum()))
    return x,var,cells,bars,qc,review


def prepare(cfg,destination):
    if cfg['study_id']!='GSE295527' or cfg['raw_location']!='X':raise ValueError('Only reviewed benchmark native RNA is enabled')
    frames=[];pins=[];first_ids=None;first_symbols=None
    for spec in cfg['benchmark_samples']:
        asset=spec['asset'];p=download(asset['url'],ROOT/'data/downloads/GSE295527_assets'/asset['filename'],asset['expected_bytes'])
        if digest(p)!=asset['sha256']:raise ValueError('Original archive checksum changed')
        x,var,source,bars,qc,review=read_capture(p,spec['platform_kind'])
        if spec['platform_kind']=='Parse' and (source['sample'].nunique()!=1 or not source['sample'].str.endswith('-'+spec['source_sample_code']).all()):
            raise ValueError('Original cell sample label conflicts with original GEO specimen')
        if any(value!=spec['source_count_review'].get(key) for key,value in review.items()):raise ValueError('Native RNA audit changed')
        qc_pin=spec['source_count_review']['original_RNA_QC']
        if digest(ROOT/qc_pin['path'])!=qc_pin['sha256']:raise ValueError('Retained original RNA QC changed')
        pd.testing.assert_frame_equal(pd.read_csv(ROOT/qc_pin['path']),qc,check_dtype=False)
        if first_ids is None:first_ids=var.index.tolist();first_symbols=var.gene_symbols.tolist()
        if first_ids!=var.index.tolist() or first_symbols!=var.gene_symbols.tolist():raise ValueError('Different native gene axes cannot be concatenated')
        # Preserve every source cell metadata column; actual project RNA QC is
        # applied once in the common pipeline after count-preserving conversion.
        obs=source.add_prefix('source_cell_').copy();obs.index=pd.Index(spec['accession']+':'+bars.astype(str),name='cell_id')
        obs['source_sample_id']=spec['accession'];obs['donor_id']=spec['donor_id'];obs['tissue']='blood';obs['is_primary_data']=True
        for key,value in spec['metadata'].items():obs[key]=value
        frames.append(obs);pins.append((spec,p))
        del x
    combined=pd.concat(frames,join='outer')
    for column in combined:
        if combined[column].dtype==object:combined[column]=pd.Categorical(combined[column])
    if Path(destination).exists():
        # An interrupted replay can leave a partial alongside the original
        # validated carrier. Preserve both and reuse the carrier only after an
        # exact count comparison with every checksum-verified source capture.
        existing=ad.read_h5ad(destination,backed='r')
        try:
            if (existing.shape!=(len(combined),len(first_ids)) or
                    existing.obs_names.tolist()!=combined.index.tolist() or
                    existing.var_names.tolist()!=first_ids):
                raise ValueError('Existing benchmark carrier axes differ')
            offset=0
            for spec,p in pins:
                x,_,_,_,_,_=read_capture(p,spec['platform_kind'])
                for start in range(0,x.shape[0],256):
                    stop=min(start+256,x.shape[0])
                    if (sparse.csr_matrix(existing.X[offset+start:offset+stop])-
                            sparse.csr_matrix(x[start:stop])).nnz:
                        raise ValueError(f'Existing benchmark carrier counts differ at {spec["accession"]} row {start}')
                offset+=x.shape[0]
        finally:
            existing.file.close()
        return dict(kind='native_benchmark_archive_RNA',assets=[s['asset'] for s,p in pins],
                    conversion='Original RNA integer counts and native gene axes unchanged; existing carrier reverified against every pinned capture'),[p for s,p in pins]
    part=Path(str(destination)+'.part')
    if part.exists():raise ValueError('Inspect previous conversion partial before retrying')
    if footprint()+2_000_000_000>=LIMIT:raise RuntimeError('Project budget exceeded')
    var=pd.DataFrame(dict(gene_symbols=first_symbols),index=first_ids)
    ad.AnnData(sparse.csr_matrix((len(combined),len(var)),dtype=np.int64),obs=combined,var=var).write_h5ad(part,compression='gzip')
    offset=nnz=0
    with h5py.File(part,'r+') as h:
        matrix=h['X']
        for key in ['data','indices','indptr']:del matrix[key]
        data=matrix.create_dataset('data',shape=(0,),maxshape=(None,),dtype='<i8',compression='gzip')
        indices=matrix.create_dataset('indices',shape=(0,),maxshape=(None,),dtype='<i4',compression='gzip')
        indptr=matrix.create_dataset('indptr',shape=(len(combined)+1,),dtype='<i8',compression='gzip');indptr[0]=0
        for spec,p in pins:
            x,_,_,_,_,review=read_capture(p,spec['platform_kind'])
            required=x.nnz*12+x.indptr.nbytes+2_000_000_000
            if footprint()+required>=LIMIT or shutil.disk_usage(ROOT).free<required:raise RuntimeError('Count conversion exceeds available disk budget')
            stop=nnz+x.nnz;data.resize((stop,));indices.resize((stop,));data[nnz:stop]=x.data;indices[nnz:stop]=x.indices
            indptr[offset+1:offset+x.shape[0]+1]=x.indptr.astype(np.int64)[1:]+nnz
            offset+=x.shape[0];nnz=stop
            print('Converted original archive',spec['accession'],x.shape,flush=True)
        if offset!=len(combined) or int(indptr[-1])!=nnz:raise ValueError('Incomplete CSR conversion')
    check=ad.read_h5ad(part,backed='r')
    try:
        if check.obs_names.tolist()!=combined.index.tolist() or check.var_names.tolist()!=first_ids:raise ValueError('Staged axes differ')
    finally:check.file.close()
    part.replace(destination)
    return dict(kind='native_benchmark_archive_RNA',assets=[s['asset'] for s,p in pins],
                conversion='Original RNA integer counts and native gene axes unchanged; Parse cell metadata count sums independently checked; no gene filtering or outer zero fill'),[p for s,p in pins]
