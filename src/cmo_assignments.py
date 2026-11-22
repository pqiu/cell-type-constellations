"""Use deposited CMO singlet calls with explicit donor mapping and threshold."""
import copy
import gzip
import io
import json
import tarfile
import numpy as np
import pandas as pd
from scripts.download_pbmc_asset import ROOT, download
from src.geo_inputs import convert_geo_tar

def select_assignments(frame, mapping, threshold):
    if not 0.5 < threshold <= 1: raise ValueError('Invalid confidence threshold')
    if not mapping.CMO_id.is_unique or not mapping.sample_id.is_unique: raise ValueError('Ambiguous CMO map')
    if not frame.Barcodes.is_unique: raise ValueError('Duplicate barcodes')
    cols=list(mapping.CMO_id)+['Multiplet','Blanks']
    probabilities=frame[cols].to_numpy(dtype=float)
    if not np.isfinite(probabilities).all() or (probabilities < -1e-12).any() or (probabilities > 1+1e-12).any(): raise ValueError('Invalid probabilities')
    if not np.allclose(probabilities.sum(axis=1),1,atol=1e-6): raise ValueError('Probabilities do not sum to one')
    called=frame.Assignment.isin(mapping.CMO_id)
    if not frame.loc[called,cols].idxmax(axis=1).equals(frame.loc[called,'Assignment']): raise ValueError('Assignment differs from probability maximum')
    expected=frame.loc[called,cols].max(axis=1).to_numpy()
    if not np.allclose(expected,frame.loc[called,'Assignment_Probability'].to_numpy(),atol=1e-6): raise ValueError('Reported confidence differs')
    selected=frame.loc[called & frame.Assignment_Probability.ge(threshold)].copy().set_index('Barcodes')
    indexed=mapping.set_index('CMO_id')
    selected['donor']=selected.Assignment.map(indexed.sample_id)
    selected['disease']=selected.Assignment.map(indexed.disease)
    if set(selected.donor)!=set(mapping.sample_id): raise ValueError('Expected donor missing')
    return selected

def prepare(cfg,destination):
    path=ROOT/'data/downloads'/(cfg['dataset_id']+'_RAW.tar')
    download(cfg['url'],path,cfg['expected_bytes'])
    derived=copy.deepcopy(cfg);overrides={};audit={}
    with tarfile.open(path) as t:
        for spec in derived['geo_samples']:
            def read(name): return pd.read_csv(io.BytesIO(gzip.decompress(t.extractfile(name).read())))
            selected=select_assignments(read(spec['confidence_member']),read(spec['cmo_map_member']),cfg['assignment_min_probability'])
            barcodes=set(gzip.decompress(t.extractfile(spec['barcodes_member']).read()).decode().splitlines())
            missing=selected.index[~selected.index.isin(barcodes)]
            if sorted(missing)!=sorted(spec['expected_missing_assignment_barcodes']): raise ValueError('RNA/CMO barcode mismatch changed')
            selected=selected.loc[selected.index.isin(barcodes)]
            member=spec['accession']+'_reviewed_assignments.csv';overrides[member]=selected.to_csv().encode()
            spec['cell_assignment']=dict(member=member,barcode_prefix='',donor_column='donor',expected_donors=cfg['expected_donors'])
            audit[spec['accession']]=len(selected)
    convert_geo_tar(path,destination,derived,metadata_overrides=overrides)
    return dict(kind='published_CMO_assignments',source=json.loads(path.with_suffix('.tar.manifest.json').read_text()),threshold=cfg['assignment_min_probability'],selected_cells=audit),[path]
