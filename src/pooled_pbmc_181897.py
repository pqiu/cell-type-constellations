"""Stream original pooled RNA while verifying every author-cell count join."""
import gzip
import hashlib
import json
import shutil
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from anndata.utils import make_index_unique
from scipy import sparse

from scripts.download_pbmc_asset import ROOT,LIMIT,footprint,download


def sha(path):
    with Path(path).open('rb') as h:return hashlib.file_digest(h,'sha256').hexdigest()


def recover(raw_path,metadata_path,destination,block_columns=131072):
    """Preserve full native RNA; accept mapping only after exact count equality.

    The source has repeated full droplet whitelists in its aggregated libraries.
    Both that structure and the candidate zero/one-based suffix conversion are
    verified; they are not accepted merely from dimensions or naming convention.
    """
    source=ad.read_h5ad(metadata_path,backed='r');obs=source.obs.copy();meta_var=source.var.copy()
    part=Path(str(destination)+'.part')
    if destination.exists() or part.exists():
        source.file.close();raise ValueError('Inspect existing RNA staging before recovery')
    try:
        if not obs.index.is_unique or not meta_var.index.is_unique:raise ValueError('Repeated metadata coordinates')
        batch=obs.batch.astype(str);bare=obs.index.str.slice(0,16)
        if not (bare+'-'+batch+'-'+batch).equals(pd.Index(obs.index)):
            # Index/Series containers differ; compare literal values explicitly.
            if not np.array_equal(np.asarray(bare+'-'+batch+'-'+batch),obs.index.to_numpy()):raise ValueError('Unreviewed author-cell suffix format')
        with h5py.File(raw_path) as raw:
            m=raw['matrix'];f=m['features'];ids=pd.Index(f['id'].asstr()[()]);symbols=f['name'].asstr()[()]
            libraries=[v.decode() if isinstance(v,bytes) else str(v) for v in raw.attrs['library_ids']]
            if libraries!=[str(i+1) for i in range(len(libraries))] or set(batch)!=set(str(i) for i in range(len(libraries))):raise ValueError('Source library/batch keys differ')
            ng,nc=map(int,m['shape'][()]);nnz=len(m['data']);per_library=nc//len(libraries)
            if nc%len(libraries) or ng!=len(ids) or not ids.is_unique or set(f['feature_type'].asstr()[()])!={'Gene Expression'} or set(f['genome'].asstr()[()])!={'GRCh38'}:raise ValueError('Unreviewed native RNA axes')
            rna=meta_var.genome.eq('GRCh38').to_numpy();projection=ids.get_indexer(meta_var.loc[rna,'gene_ids'])
            # Source H5AD uniquifies duplicated symbols with AnnData's '-' join.
            # Verify the full native-axis transformation and authoritative IDs;
            # do not strip arbitrary suffixes or replace native gene symbols.
            unique_symbols=make_index_unique(pd.Index(symbols),join='-')
            if (projection<0).any() or not np.array_equal(unique_symbols[projection],meta_var.index[rna].to_numpy()):raise ValueError('Processed/native RNA identities differ')
            renamed_symbols=int(np.count_nonzero(symbols[projection]!=meta_var.index[rna].to_numpy()))
            if not np.array_equal(np.concatenate([np.flatnonzero(batch.eq(str(i))) for i in range(len(libraries))]),np.arange(len(obs))):raise ValueError('Metadata library order differs; explicit permutation required')
            needed=nnz*12+2_000_000_000
            if footprint()+needed>=LIMIT or shutil.disk_usage(ROOT).free<needed:raise RuntimeError('Insufficient project budget for native RNA recovery')
            var=pd.DataFrame(dict(gene_symbols=symbols),index=ids)
            ad.AnnData(sparse.csr_matrix((len(obs),ng),dtype=np.int32),obs=obs,var=var).write_h5ad(part,compression='gzip')
            offset=written=previous_pointer=total=zero_entries=0;qc_parts=[];library_reviews=[]
            with h5py.File(part,'r+') as out:
                mat=out['X']
                for key in ['data','indices','indptr']:del mat[key]
                values_out=mat.create_dataset('data',shape=(0,),maxshape=(None,),dtype='<i4',compression='gzip')
                indices_out=mat.create_dataset('indices',shape=(0,),maxshape=(None,),dtype='<i4',compression='gzip')
                pointers_out=mat.create_dataset('indptr',shape=(len(obs)+1,),dtype='<i8',compression='gzip');pointers_out[0]=0
                for lib_index,lib in enumerate(libraries):
                    wanted_rows=np.flatnonzero(batch.eq(str(lib_index)));wanted=np.asarray(bare[wanted_rows],dtype='S16');cursor=0;last_base=None
                    if len(wanted)>1 and not (wanted[1:]>wanted[:-1]).all():raise ValueError('Metadata barcodes are not unique and ordered within a library')
                    for local in range(0,per_library,block_columns):
                        stop_local=min(local+block_columns,per_library);start=lib_index*per_library+local;stop=lib_index*per_library+stop_local
                        bars=m['barcodes'][start:stop];base=bars.astype('S16')
                        if lib_index==0:
                            chars=base.view('u1').reshape(-1,16)
                            if not np.isin(chars,[65,67,71,84]).all():raise ValueError('Invalid original droplet bases')
                        elif not np.array_equal(base,m['barcodes'][local:stop_local].astype('S16')):raise ValueError('Native repeated-whitelist structure changed')
                        if not np.array_equal(bars,np.char.add(base,('-'+lib).encode())) or (len(base)>1 and not (base[1:]>base[:-1]).all()) or (last_base is not None and base[0]<=last_base):raise ValueError('Repeated or malformed native barcodes')
                        last_base=base[-1]
                        ptr=m['indptr'][start:stop+1];data=m['data'][int(ptr[0]):int(ptr[-1])];ix=m['indices'][int(ptr[0]):int(ptr[-1])]
                        if ptr[0]!=previous_pointer or (np.diff(ptr)<0).any() or ptr[-1]>nnz:raise ValueError('Invalid full native CSC pointers')
                        previous_pointer=int(ptr[-1])
                        if not np.isfinite(data).all() or (data<0).any() or not np.equal(data,np.floor(data)).all() or (ix<0).any() or (ix>=ng).any():raise ValueError('Invalid full native values or row coordinates')
                        bad=np.flatnonzero(ix[1:]<=ix[:-1])+1
                        if len(bad) and not np.isin(bad,ptr-ptr[0]).all():raise ValueError('Repeated or unsorted within-cell native gene coordinates')
                        total+=int(data.sum(dtype=np.int64));zero_entries+=int(np.count_nonzero(data==0))
                        right=int(np.searchsorted(wanted,base[-1],side='right'));selected=np.searchsorted(base,wanted[cursor:right])
                        if not np.array_equal(base[selected],wanted[cursor:right]):raise ValueError('Author called-cell barcode missing in raw RNA')
                        if len(selected):
                            rows=wanted_rows[cursor:right]
                            if not np.array_equal(rows,np.arange(offset,offset+len(rows))):raise ValueError('Author-cell row order differs from verified native order')
                            x=sparse.csc_matrix((data,ix,ptr-ptr[0]),shape=(ng,stop-start))[:,selected].T.tocsr();x.eliminate_zeros()
                            deposited=source.X[rows].tocsr()[:,rna].tocsr();deposited.sum_duplicates();deposited.sort_indices();deposited.eliminate_zeros()
                            projected=x[:,projection].tocsr();projected.sort_indices()
                            for key in ['indptr','indices','data']:
                                if not np.array_equal(getattr(projected,key),getattr(deposited,key)):raise ValueError('Candidate barcode mapping fails full original RNA count equality')
                            end=written+x.nnz;values_out.resize((end,));indices_out.resize((end,));values_out[written:end]=x.data;indices_out[written:end]=x.indices
                            pointers_out[offset+1:offset+len(rows)+1]=x.indptr.astype(np.int64)[1:]+written
                            libraries_total=np.asarray(x.sum(axis=1,dtype=np.int64)).ravel();filtered_total=np.asarray(deposited.sum(axis=1,dtype=np.int64)).ravel()
                            mito=np.asarray(x[:,np.array([s.startswith('MT-') for s in symbols])].sum(axis=1,dtype=np.int64)).ravel()
                            qc_parts.append(pd.DataFrame(dict(source_cell_id=obs.index[rows],source_raw_barcode=[v.decode() for v in bars[selected]],source_raw_column=start+selected,RNA_total=libraries_total,RNA_detected_genes=np.diff(x.indptr),mito_RNA_total=mito,deposited_filtered_RNA_total=filtered_total)))
                            offset+=len(rows);written=end
                        cursor=right
                    if cursor!=len(wanted):raise ValueError('Not every author cell recovered')
                    library_reviews.append(dict(original_library=lib,author_batch=str(lib_index),raw_droplet_whitelist_rows=per_library,exact_count_matched_author_cells=cursor))
                    print('Verified original library',lib,cursor,'author cells; full native values/coordinates checked',flush=True)
                if previous_pointer!=nnz or offset!=len(obs) or pointers_out[-1]!=written:raise ValueError('Incomplete original matrix stream')
            qc=pd.concat(qc_parts,ignore_index=True)
            if qc.source_cell_id.tolist()!=obs.index.tolist() or qc.source_raw_barcode.duplicated().any():raise ValueError('Recovered source cell identities differ')
            part.replace(destination)
            return dict(full_original_count_and_coordinate_audit=True,source_raw_droplet_rows=nc,source_stored_values=nnz,source_zero_entries=zero_entries,full_source_RNA_count_sum=total,called_cells=len(obs),native_RNA_genes=ng,matched_deposited_RNA_genes=len(projection),native_mito_genes=sum(s.startswith('MT-') for s in symbols),all_native_barcodes_unique=True,all_author_cells_match_exact_native_RNA=True,source_uniquified_gene_symbols=renamed_symbols,source_symbol_conversion='Exact AnnData make_index_unique(join=-) on full native symbols matches source H5AD by authoritative Ensembl ID. Full original symbols retained as geometry annotations.',native_RNA_gene_mapping='Original Ensembl IDs and symbols retained;20303deposited RNA genes exactly project from33538native RNA features.',barcode_mapping='Original16base-(one-based library) equals author16base-(zero-based batch)-(same batch); all original barcodes and every available RNA count verify the mapping.',libraries=library_reviews,derived_stage_sha256=sha(destination)),qc
    finally:source.file.close()


def expand_metadata(path):
    expanded=path.with_suffix('');manifest=expanded.with_suffix('.h5ad.manifest.json')
    if expanded.exists():
        review=json.loads(manifest.read_text())
        if review['compressed_source']['sha256']!=sha(path) or review['sha256']!=sha(expanded):raise ValueError('Expanded metadata bytes changed')
        return expanded
    part=Path(str(expanded)+'.part')
    if part.exists():raise ValueError('Inspect partial metadata expansion')
    allowance=min(LIMIT-footprint(),shutil.disk_usage(ROOT).free)-2_000_000_000;size=0
    with gzip.open(path,'rb') as source,part.open('wb') as out:
        while block:=source.read(8*1024*1024):
            size+=len(block)
            if size>allowance:raise RuntimeError('Metadata expansion exceeds project budget')
            out.write(block)
    part.replace(expanded);manifest.write_text(json.dumps(dict(kind='unchanged_gzip_expansion',compressed_source=dict(path=path.relative_to(ROOT).as_posix(),sha256=sha(path)),bytes=size,sha256=sha(expanded)),indent=2)+'\n')
    return expanded


def prepare(cfg,destination):
    if cfg['study_id']!='GSE181897':raise ValueError('Explicit pooled source review required')
    paths={}
    for role,asset in cfg['source_assets'].items():
        if Path(asset['filename']).name!=asset['filename']:raise ValueError('Unsafe original source name')
        p=download(asset['url'],ROOT/'data/downloads'/asset['filename'],asset['expected_bytes'])
        if sha(p)!=asset['sha256']:raise ValueError('Original source checksum changed')
        paths[role]=p
    metadata=expand_metadata(paths['metadata']);pin=cfg['original_RNA_recovery_review'];audit_path=ROOT/pin['path']
    if sha(audit_path)!=pin['sha256']:raise ValueError('Original RNA recovery review changed')
    audit=json.loads(audit_path.read_text())
    if destination.exists():
        if sha(destination)!=audit['derived_stage_sha256']:raise ValueError('Existing native RNA stage differs from reviewed bytes')
        review=audit
    else:
        review,qc=recover(paths['RNA'],metadata,destination)
        for key in ['full_source_RNA_count_sum','called_cells','source_stored_values','all_author_cells_match_exact_native_RNA','native_RNA_genes']:
            if review[key]!=audit[key]:raise ValueError('Original RNA source recovery changed')
        pin=audit['retained_native_QC']
        if sha(ROOT/pin['path'])!=pin['sha256']:raise ValueError('Original RNA QC pin changed')
        pd.testing.assert_frame_equal(pd.read_csv(ROOT/pin['path']),qc,check_dtype=False)
    return dict(kind='fully_verified_original_pooled_RNA',recovery=review,source_assets=cfg['source_assets']),list(paths.values())+[metadata]
