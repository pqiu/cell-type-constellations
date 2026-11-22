"""Use deposited Vireo calls and explicit GEO genetic-cluster crosswalks."""
import copy
import hashlib
import json
import re
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scripts.download_pbmc_asset import ROOT
from src.geo_assets import prepare_assets
from src.geo_inputs import convert_geo_tar


def select(frame, mapping, threshold=0.90):
    if not 0.5<=threshold<=1:raise ValueError('Invalid assignment threshold')
    if not frame.cell.is_unique or frame.cell.isna().any():raise ValueError('Duplicate/missing Vireo cell identity')
    if set(frame.donor_id)-set(mapping)-{'doublet','unassigned'}:raise ValueError('Unknown Vireo call')
    p=frame[['prob_max','prob_doublet']].to_numpy(dtype=float)
    if not np.isfinite(p).all() or (p<0).any() or (p>1).any():raise ValueError('Invalid Vireo probability')
    keep=frame.donor_id.isin(mapping)&frame.prob_max.ge(threshold)&frame.prob_doublet.le(1-threshold+1e-12)&frame.n_vars.gt(0)
    result=frame.loc[keep].copy().set_index('cell')
    if not result.best_singlet.eq(result.donor_id).all():raise ValueError('Vireo call disagrees with best singlet')
    if set(result.donor_id)!=set(mapping):raise ValueError('Source contributor lacks retained assignments')
    result['specimen_id']=result.donor_id.map(mapping)
    if not result.specimen_id.str.fullmatch(r'L\d+_T[14]').all():raise ValueError('Unreviewed source specimen label')
    result['biological_donor']=result.specimen_id.str.replace(r'_T[14]$','',regex=True)
    result['source_visit']=result.specimen_id.str.extract(r'_(T[14])$',expand=False)
    return result


def audit_h5(path, selected, seen, legacy_genome=None):
    with h5py.File(path,'r') as f:
        if legacy_genome is not None:
            if list(f)!=[legacy_genome]:raise ValueError('Legacy genome group differs from reviewed source')
            m=f[legacy_genome];ids=m['genes'][()];symbols=m['gene_names'][()]
            types=np.full(len(ids),b'Gene Expression')
        else:
            m=f['matrix'];ft=m['features'];ids=ft['id'][()];symbols=ft['name'][()];types=ft['feature_type'][()]
        barcodes=[x.decode() for x in m['barcodes'][()]]
        if len(barcodes)!=len(set(barcodes)) or not selected.index.isin(barcodes).all():raise ValueError('Vireo/RNA cell identity mismatch')
        v=m['data'][()]
        if not np.isfinite(v).all() or (v<0).any() or not np.equal(v,np.floor(v)).all():raise ValueError('Invalid original UMI counts')
        x=sparse.csc_matrix((v,m['indices'][()],m['indptr'][()]),shape=tuple(m['shape'][()])).T.tocsr()
        x.check_format(full_check=True)
        if x.shape!=(len(barcodes),len(ids)):raise ValueError('Source matrix identity dimensions mismatch')
        rna=types==b'Gene Expression'
        if not rna.any() or len(ids[rna])!=len(set(ids[rna])):raise ValueError('Invalid original RNA coordinates')
        x=x[:,rna];x.sum_duplicates();x.sort_indices();x.eliminate_zeros()
        gene_sha=hashlib.sha256(json.dumps([[v.decode() for v in ids[rna]],[v.decode() for v in symbols[rna]]],separators=(',',':')).encode()).hexdigest()
        lookup={b:i for i,b in enumerate(barcodes)}
        for barcode in selected.index:
            i=lookup[barcode];start,end=x.indptr[i:i+2]
            if end-start<20:continue
            fingerprint=hashlib.sha256(x.indices[start:end].astype('<i8').tobytes()+x.data[start:end].astype('<i8').tobytes()).hexdigest()
            key=(barcode,gene_sha,fingerprint)
            if key in seen:raise ValueError('Exact source cell/count duplication across libraries requires review')
            seen.add(key)
        return dict(shape=list(x.shape),nonzeros=x.nnz,integer_nonnegative_counts=True,gene_coordinate_sha256=gene_sha)


def prepare(cfg,destination):
    archive,manifest,paths=prepare_assets(cfg)
    derived=copy.deepcopy(cfg);overrides={};reviews={};seen=set();genes=set();specimens=set();donors=set()
    for spec in derived['geo_samples']:
        path=ROOT/f'data/downloads/{cfg["study_id"]}_assets'/spec['vireo_member']
        selected=select(pd.read_csv(path,sep='\t'),spec['vireo_mapping'],cfg['assignment_min_probability'])
        source=ROOT/f'data/downloads/{cfg["study_id"]}_assets'/spec['matrix_h5_member']
        review=audit_h5(source,selected,seen);genes.add(review['gene_coordinate_sha256'])
        specimens.update(selected.specimen_id);donors.update(selected.biological_donor)
        name=spec['accession']+'_reviewed_assignments.csv';overrides[name]=selected.to_csv().encode()
        spec['cell_assignment']=dict(member=name,barcode_prefix='',donor_column='biological_donor',expected_donors=sorted(selected.biological_donor.unique()))
        reviews[spec['accession']]=dict(**review,selected_cells=len(selected),specimens=sorted(selected.specimen_id.unique()))
    if len(genes)!=1 or sorted(specimens)!=cfg['expected_specimens'] or sorted(donors)!=cfg['expected_donors']:
        raise ValueError('Changed gene coverage or biological identity universe')
    convert_geo_tar(archive,destination,derived,metadata_overrides=overrides)
    return dict(kind='source_Vireo_GEO_crosswalk',assets=manifest,assignment_min_probability=cfg['assignment_min_probability'],
        source_libraries=reviews,uncertainty='Vireo model-based confidence, not calibrated ground truth; doublet/unassigned calls excluded'),[archive]+paths
