"""Load checksum-pinned author count members with retained cell assignments."""
import gzip
import hashlib
import json
import tarfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
from threadpoolctl import threadpool_limits

from scripts.download_pbmc_asset import ROOT, LIMIT, download, footprint


def sha256(path):
    with Path(path).open('rb') as h:
        return hashlib.file_digest(h,'sha256').hexdigest()


def download_path(relative):
    path=(ROOT/relative).resolve()
    if not path.is_relative_to((ROOT/'data/downloads').resolve()):
        raise ValueError('Count carriers must remain inside project downloads')
    return path


def prepare(cfg,destination):
    destination=Path(destination).resolve()
    if not destination.is_relative_to((ROOT/'data/downloads').resolve()):
        raise ValueError('Converted count output must remain inside project downloads')
    archive_pin=cfg['source_archive']
    archive=download_path(archive_pin['path'])
    download(archive_pin['url'],archive,archive_pin['bytes'])
    if sha256(archive)!=archive_pin['sha256']:
        raise ValueError('Original author archive checksum changed')
    expected={}
    for capture in cfg['author_captures']:
        for asset in capture['source_assets']:
            if asset['member'] in expected:
                raise ValueError('Duplicate source archive member')
            download_path(asset['path'])
            expected[asset['member']]=asset
    if not expected:
        raise ValueError('Reviewed source count members required')
    missing=[]
    for member,pin in expected.items():
        path=download_path(pin['path'])
        if path.exists():
            if path.stat().st_size!=pin['bytes'] or sha256(path)!=pin['sha256']:
                raise ValueError('Retained original count member changed')
        else: missing.append(member)
    if footprint()+sum(expected[n]['bytes'] for n in missing)+2_000_000_000>=LIMIT:
        raise RuntimeError('Author count extraction exceeds project budget')
    if missing:
        found=set()
        with tarfile.open(archive,'r|*') as tar:
            for member in tar:
                if member.name not in expected: continue
                if member.name in found or not member.isfile() or member.size!=expected[member.name]['bytes']:
                    raise ValueError('Invalid or duplicate original count member')
                found.add(member.name)
                if member.name not in missing: continue
                pin=expected[member.name]; path=download_path(pin['path']); path.parent.mkdir(parents=True,exist_ok=True)
                part=path.with_suffix(path.suffix+'.extracting'); h=hashlib.sha256()
                with tar.extractfile(member) as source,part.open('wb') as out:
                    while block:=source.read(8*1024*1024): out.write(block); h.update(block)
                if part.stat().st_size!=pin['bytes'] or h.hexdigest()!=pin['sha256']:
                    raise ValueError('Extracted original count member checksum changed')
                part.replace(path)
        if found!=set(expected): raise ValueError('Missing original archive members')
    pieces=[]; tables=[]; first_ids=None; first_symbols=None; evidence=[]
    for capture in cfg['author_captures']:
        assets={Path(p['path']).name:download_path(p['path']) for p in capture['source_assets']}
        name=capture['capture']
        genes=pd.read_csv(assets[name+'_features.tsv.gz'],sep='\t',header=None,dtype=str,keep_default_na=False)
        barcodes=pd.read_csv(assets[name+'_barcodes.tsv.gz'],sep='\t',header=None,dtype=str,keep_default_na=False)[0]
        if genes.shape[1]!=3 or not genes[0].is_unique or not barcodes.is_unique:
            raise ValueError('Invalid native count coordinates')
        if genes[0].str.strip().eq('').any() or genes[1].str.strip().eq('').any() or barcodes.str.strip().eq('').any():
            raise ValueError('Empty native gene or barcode coordinate')
        if not set(genes[2]).issubset({'Gene Expression','Antibody Capture'}):
            raise ValueError('Unreviewed source feature modality')
        with threadpool_limits(8),gzip.open(assets[name+'_matrix.mtx.gz'],'rb') as h:
            x=mmread(h,spmatrix=True).tocsc()
        x.check_format(full_check=True)
        if x.shape!=(len(genes),len(barcodes)) or not np.isfinite(x.data).all() or (x.data<0).any() or not np.equal(x.data,np.floor(x.data)).all():
            raise ValueError('Invalid original nonnegative integer counts')
        x.sum_duplicates(); x.sort_indices(); x.eliminate_zeros()
        rna=genes[2].eq('Gene Expression').to_numpy()
        if not rna.any(): raise ValueError('No native RNA features')
        native=genes.loc[rna].reset_index(drop=True)
        ids=native[0].to_numpy(); symbols=native[1].to_numpy()
        if first_ids is None: first_ids=ids; first_symbols=symbols
        if not np.array_equal(ids,first_ids) or not np.array_equal(symbols,first_symbols):
            raise ValueError('Differing native gene axes need explicit availability handling')
        raw=x[rna].T.tocsr(); libraries=np.asarray(raw.sum(axis=1)).ravel()
        if int(raw.sum())!=capture['source_total_RNA_counts']:
            raise ValueError('Original RNA-only total changed')
        pin=capture['assignment']; assignment=(ROOT/pin['path']).resolve()
        if not assignment.is_relative_to((ROOT/'references/datasets').resolve()) or sha256(assignment)!=pin['sha256']:
            raise ValueError('Retained cell assignment checksum or location changed')
        frame=pd.read_csv(assignment,index_col=0,dtype=str,keep_default_na=False)
        if not frame.index.is_unique or not frame.index.isin(barcodes).all() or not frame.source_capture.eq(name).all():
            raise ValueError('Assignment/capture/native-barcode mismatch')
        keep=barcodes.isin(frame.index).to_numpy(); selected_barcodes=barcodes.loc[keep]
        frame=frame.reindex(selected_barcodes)
        for key in ['source_person_code','source_specimen_id','source_visit_code']:
            if frame[key].isna().any() or frame[key].str.strip().eq('').any():
                raise ValueError('Missing original donor/specimen/visit assignment')
        if not np.array_equal(pd.to_numeric(frame.RNA_total).to_numpy(),libraries[keep]):
            raise ValueError('Assigned original RNA library totals changed')
        selected=raw[keep].tocsr()
        if not np.array_equal(pd.to_numeric(frame.RNA_detected_genes).to_numpy(),np.diff(selected.indptr)):
            raise ValueError('Assigned original detected RNA genes changed')
        obs=frame.copy(); obs.index=pd.Index(name+':'+selected_barcodes)
        obs.index.name='cell_id'
        obs['donor_id']=frame.source_person_code.to_numpy(); obs['source_sample_id']=frame.source_specimen_id.to_numpy()
        obs['tissue']='blood'; obs['source_specimen']='PBMC'; obs['is_primary_data']=True
        for key,value in cfg.get('cell_metadata_constants',{}).items():
            if key in obs: raise ValueError('Refuse original assignment metadata overwrite')
            obs[key]=value
        pieces.append(selected); tables.append(obs)
        evidence.append(dict(source_capture=name,original_barcode_rows=len(barcodes),assigned_cells=len(obs),
            source_RNA_total=int(raw.sum()),assignment=pin))
    counts=sparse.vstack(pieces,format='csr'); obs=pd.concat(tables,axis=0)
    if not obs.index.is_unique: raise ValueError('Duplicated cross-capture cell identity')
    for column in obs:
        if obs[column].dtype==object: obs[column]=pd.Categorical(obs[column])
    var=pd.DataFrame({'gene_symbols':first_symbols,'feature_name':first_symbols},index=first_ids)
    alias_plan=cfg.get('annotation_gene_symbol_aliases')
    if alias_plan is not None:
        aliases=alias_plan['mapping'];pin=alias_plan['evidence']
        proof_path=(ROOT/pin['path']).resolve()
        if not proof_path.is_relative_to((ROOT/'references/datasets').resolve()) or sha256(proof_path)!=pin['sha256']:
            raise ValueError('Annotation alias evidence checksum or location changed')
        proof=json.loads(proof_path.read_text(encoding='utf-8'))
        if proof.get('validated') is not True or proof.get(alias_plan['evidence_mapping_key'])!=aliases:
            raise ValueError('Annotation aliases disagree with validated source evidence')
        if not isinstance(aliases,dict) or not aliases or not set(aliases).issubset(set(first_ids)):
            raise ValueError('Annotation aliases require exact present source gene IDs')
        if any(not isinstance(v,str) or not v.strip() for v in aliases.values()) or len(set(aliases.values()))!=len(aliases):
            raise ValueError('Annotation aliases require nonempty one-to-one target symbols')
        unmapped={symbol for gene,symbol in zip(first_ids,first_symbols) if gene not in aliases}
        if set(aliases.values()) & unmapped:
            raise ValueError('Annotation alias collides with an unmapped source gene')
        # Raw IDs, original symbols and all expression coordinates/values stay intact.
        var['annotation_gene_symbols']=[aliases.get(gene,symbol) for gene,symbol in zip(first_ids,first_symbols)]
    estimated=2*(counts.data.nbytes+counts.indices.nbytes+counts.indptr.nbytes)+2*int(obs.memory_usage(deep=True).sum())+100_000_000
    if footprint()+estimated>=LIMIT: raise RuntimeError('Converted author counts exceed project budget')
    obj=ad.AnnData(X=counts,obs=obs,var=var)
    obj.uns['input_note']='Pinned original RNA counts in X; antibody counts excluded; retained derived cell assignments and uncertainty preserved.'
    obj.uns['source_preconversion_audit_json']=json.dumps(evidence)
    obj.write_h5ad(destination,compression='gzip')
    manifest=dict(kind='pinned_original_author_count_members',source_archive=archive_pin,
        native_count_members=list(expected.values()),capture_assignments=evidence)
    if alias_plan is not None:
        manifest['annotation_gene_symbol_aliases']=alias_plan
    return manifest,[archive]+[download_path(p['path']) for p in expected.values()]
