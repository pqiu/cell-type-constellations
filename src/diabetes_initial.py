"""Combine final experiment cells with explicit source feature availability."""
import hashlib
import json
from pathlib import Path
import zipfile
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from src.rds_seurat import read_count_slots
from src.rdata_counts import sparse_counts
from scripts.download_pbmc_asset import ROOT, LIMIT, footprint, download

def json_value(value):
    if value is None or value is pd.NA or isinstance(value,np.ma.core.MaskedConstant):
        return None
    if isinstance(value,np.ndarray):
        return json_value(value.tolist())
    if isinstance(value,np.generic):
        return json_value(value.item())
    if isinstance(value,(list,tuple)):
        return [json_value(item) for item in value]
    if isinstance(value,(float,int,bool,str)):
        return None if isinstance(value,float) and not np.isfinite(value) else value
    raise ValueError('Unsupported source metadata value: '+type(value).__name__)

def prepare(cfg,destination):
    asset = cfg['source_asset']
    name = asset['filename']
    if Path(name).name != name:
        raise ValueError('Unsafe original archive name')
    archive = download(asset['url'],ROOT/'data/downloads'/name,asset['expected_bytes'])
    manifest = json.loads(archive.with_suffix(archive.suffix+'.manifest.json').read_text())
    folder = ROOT/'data/downloads/GSE309970_assets'
    folder.mkdir(exist_ok=True)
    parts=[]; paths=[]; source_reviews=[]
    with zipfile.ZipFile(archive) as handle:
        for library,review in enumerate(cfg['initial_member_reviews']):
            member = review['source_member']
            info = handle.getinfo(member)
            if info.file_size > 500_000_000 or footprint()+info.file_size+2_000_000_000>=LIMIT:
                raise RuntimeError('Source expansion exceeds reviewed budget')
            path = folder/Path(member).name
            content = handle.read(info)
            sha = hashlib.sha256(content).hexdigest()
            if sha != review['source_sha256']:
                raise ValueError('Pinned source member changed')
            if path.exists() and hashlib.sha256(path.read_bytes()).hexdigest()!=sha:
                raise ValueError('Existing extraction differs')
            if not path.exists():
                path.write_bytes(content)
            del content
            raw,metadata = read_count_slots(path,'RNA')
            counts,genes,cells = sparse_counts(raw)
            if list(counts.shape)!=review['shape'] or hashlib.sha256('\n'.join(genes).encode()).hexdigest()!=review['gene_sha256']:
                raise ValueError('Pinned RNA dimensions/coordinates changed')
            if not metadata.index.is_unique or set(metadata.index)!=set(cells):
                raise ValueError('Metadata/count barcode mismatch')
            obs=metadata.loc[cells].copy()
            totals=np.asarray(counts.sum(axis=1)).ravel()
            if not np.array_equal(totals,obs.nCount_RNA.to_numpy()):
                raise ValueError('RNA source library total differs')
            for column in obs:
                if isinstance(obs[column].dtype,pd.CategoricalDtype):
                    obs[column]=obs[column].astype(object)
                if obs[column].dtype==object:
                    def convert(value):
                        converted=json_value(value)
                        if isinstance(converted,list):
                            return json.dumps(converted,separators=(',',':'))
                        return 'not_reported' if converted is None else str(converted)
                    obs[column]=obs[column].map(convert)
            for field in ['Patient_ID','Condition','Sample_ID']:
                if field not in obs or obs[field].astype(str).isin(['not_reported','nan','']).any():
                    raise ValueError('Unresolved source identity')
                obs[field]=obs[field].astype(str)
            obs['source_barcode']=cells.to_numpy()
            obs.index=Path(member).stem+':'+cells
            obs['source_library_index']=library
            obs['source_library_id']=Path(member).stem
            obs['donor_id']=obs.Patient_ID
            obs['tissue']='blood';obs['is_primary_data']=True
            obs['source_specimen']='PBMC-derived sorted CD4/CD8 T cells'
            obs['source_sort_enrichment']='Final experiment:1000naive and5000antigen-experiencedcells perCD4/CD8sort; notwholePBMCcomposition'
            parts.append((counts,genes,obs));paths.append(path)
            source_reviews.append(dict(source_member=member,sha256=sha,shape=list(counts.shape),gene_sha256=review['gene_sha256']))
    if len(parts)!=16:
        raise ValueError('Exactly16finalexperimentwells required')
    union=pd.Index(sorted({gene for counts,genes,obs in parts for gene in genes}))
    availability=np.zeros((len(parts),len(union)),dtype=bool)
    matrices=[]
    for library,(counts,genes,obs) in enumerate(parts):
        columns=union.get_indexer(genes)
        if (columns<0).any():
            raise ValueError('Invalid union mapping')
        availability[library,columns]=True
        coo=counts.tocoo()
        matrices.append(sparse.csr_matrix((coo.data,(coo.row,columns[coo.col])),shape=(len(obs),len(union))))
    combined=sparse.vstack(matrices,format='csr')
    if footprint()+2*(combined.data.nbytes+combined.indices.nbytes+combined.indptr.nbytes)+2_000_000_000>=LIMIT:
        raise RuntimeError('Combined conversion exceeds disk budget')
    obs=pd.concat([part[2] for part in parts])
    if not obs.index.is_unique:
        raise ValueError('Duplicate source-library cell identities')
    keys=obs.Patient_ID+'|'+obs.Condition
    if keys.nunique()!=80 or obs.donor_id.nunique()!=48:
        raise ValueError('Expected80sourcevisitgroups48donors')
    for field in ['Age','Sex']:
        if obs.groupby(keys)[field].nunique(dropna=False).max()>1:
            raise ValueError('Conflicting clinical metadata within source visit')
    var=pd.DataFrame({'gene_symbols':union.to_numpy()},index=union)
    a=ad.AnnData(X=combined,obs=obs,var=var)
    a.uns['source_gene_availability']=availability
    a.uns['source_gene_ids']=union.to_numpy()
    a.uns['source_library_names']=[Path(review['source_member']).stem for review in cfg['initial_member_reviews']]
    a.write_h5ad(destination,compression='gzip')
    return dict(kind='originalfinalexperimentRNAwithsourceavailability',archive=manifest,source_members=source_reviews,
        conversion='Union storage uses zero placeholders only with explicitavailability; incomplete all-cell gene statistics must be NaN. No expression normalization duringconversion.'),[archive]+paths
