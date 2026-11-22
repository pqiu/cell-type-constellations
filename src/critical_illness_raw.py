"""Recover original PBMC counts while retaining published called-cell metadata."""
import json
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from src.cellranger_called_cells import read_selected
from scripts.download_pbmc_asset import ROOT,LIMIT,footprint,download

def prepare(cfg,destination):
    asset=cfg['metadata_source_asset']
    if Path(asset['filename']).name!=asset['filename']:
        raise ValueError('Unsafe called-cell source name')
    source=download(asset['url'],ROOT/'data/downloads'/asset['filename'],asset['expected_bytes'])
    source_manifest=json.loads(source.with_suffix(source.suffix+'.manifest.json').read_text())
    a=ad.read_h5ad(source,backed='r')
    metadata=a.obs.copy();a.file.close()
    libraries={row['library_id'] for row in cfg['raw_matrix_sources']}
    if len(libraries)!=47 or set(metadata['orig.ident'])!=libraries or not metadata.index.is_unique:
        raise ValueError('Expected exact47called-cell library identities')
    parts=[];paths=[source];manifests=[];reference=None;var=None
    for row in cfg['raw_matrix_sources']:
        library=row['library_id'];obs=metadata.loc[metadata['orig.ident'].eq(library)].copy()
        prefix=library+'-'
        if not obs.index.str.startswith(prefix).all() or set(obs.SCARAB_ID)!={row['subject_id']}:
            raise ValueError('Called-cell and participant identifiers disagree')
        url=row['url'];name=url.rsplit('/',1)[1]
        path=download(url,ROOT/'data/downloads/GSE290679_assets'/name)
        counts,genes,symbols,review=read_selected(path,obs.index.str[len(prefix):],library)
        if reference is not None and not genes.equals(reference):
            raise ValueError('Different original RNA universes require availability handling')
        reference=genes
        if var is None:
            var=pd.DataFrame({'gene_symbols':symbols},index=genes)
        elif not np.array_equal(var.gene_symbols.to_numpy(),symbols):
            raise ValueError('Original gene symbols differ')
        original_totals=np.asarray(counts.sum(axis=1)).ravel()
        if (original_totals<obs.nCount_RNA.to_numpy()).any():
            raise ValueError('Original count totals below corrected source totals')
        obs=obs.rename(columns={column:'source_corrected_'+column for column in ['nCount_RNA','nFeature_RNA','percent.mt','percent.ribosomal','percent.heme']})
        for column in obs:
            if isinstance(obs[column].dtype,pd.CategoricalDtype):
                obs[column]=obs[column].astype(object)
            if obs[column].dtype==object:
                obs[column]=obs[column].fillna('not_reported').astype(str)
        fields=row['clinical_source']
        obs['source_sample_id']=row['accession'];obs['donor_id']=row['subject_id']
        obs['source_barcode']=obs.index.str[len(prefix):]
        obs['source_library_id']=library;obs['source_patient_id']=row['subject_id']
        obs['tissue']='blood';obs['is_primary_data']=True;obs['source_specimen']='PBMC'
        obs['source_group_geo']=fields['group'];obs['source_age_geo']=fields['age']
        obs['source_sex_geo']=fields['Sex'];obs['source_icu_geo']=fields['location']
        obs['source_age_topcoded']=float(fields['age'])>=90
        obs['author_cell_type']=obs.celltype_fine.astype(str)
        obs['disease']=obs.Group.map({'CI-Sep':'sepsis critical illness','CI-NS':'nonseptic critical illness','NHC':'non-acutely ill control'})
        if obs.disease.isna().any():
            raise ValueError('Unknown source clinical group')
        obs['source_collection_protocol']='ICU48-72hoursafteradmission; controlscommunity/bedside-surrogates, exactindividualcollectiontimeunavailable'
        parts.append((counts,obs));paths.append(path)
        manifests.append(dict(source=json.loads(path.with_suffix(path.suffix+'.manifest.json').read_text()),recovery=review,library_id=library,subject_id=row['subject_id']))
    counts=sparse.vstack([part[0] for part in parts],format='csr')
    obs=pd.concat([part[1] for part in parts])
    if len(obs)!=644147 or obs.source_sample_id.nunique()!=47 or not obs.index.is_unique:
        raise ValueError('Unexpected original called-cell coverage')
    if footprint()+2*(counts.data.nbytes+counts.indices.nbytes+counts.indptr.nbytes)+2_000_000_000>=LIMIT:
        raise RuntimeError('Original count conversion exceeds disk budget')
    ad.AnnData(X=counts,obs=obs,var=var).write_h5ad(destination,compression='gzip')
    return dict(kind='originalCellRangerRNAonpublishedcalledcells',called_cell_metadata_source=source_manifest,
        original_matrix_sources=manifests,
        conversion='CombinedsourceX isCellBendercorrected andunusedforexpression. OriginalRNAcounts selectedbyexactcalledbarcodes; noambientcorrection, embedding, ornormalization applied.'),paths
