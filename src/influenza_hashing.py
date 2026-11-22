"""Replay the reviewed influenza hashing audit before RNA count conversion."""
import hashlib
import importlib.metadata
import json
from pathlib import Path
import pandas as pd
from scripts.download_pbmc_asset import ROOT, download
from scripts.audit_influenza_hashing import audit
from src.geo_mtx_table import convert


def prepare(cfg, destination):
    if cfg['study_id'] != 'GSE309140' or importlib.metadata.version('scanpy') != cfg['hashing_scanpy_version']:
        raise ValueError('Reviewed study and hashing software required')
    crosswalk = ROOT/'references/datasets/GSE309140_Reaction_HTO_table.csv.gz'
    if hashlib.sha256(crosswalk.read_bytes()).hexdigest() != cfg['hashing_crosswalk_sha256']:
        raise ValueError('Source hashing crosswalk changed')
    paths = {}; manifests = {}
    for role, asset in cfg['matrix_assets'].items():
        if role not in ('matrix', 'features', 'barcodes') or Path(asset['filename']).name != asset['filename']:
            raise ValueError('Unexpected source asset')
        paths[role] = download(asset['url'], ROOT/'data/downloads'/asset['filename'], asset['expected_bytes'])
        manifests[role] = json.loads(paths[role].with_suffix(paths[role].suffix+'.manifest.json').read_text())
    assignments, report = audit()
    if set(assignments.source_sample_id) != {r['source_sample_id'] for r in cfg['specimen_metadata']}:
        raise ValueError('Assigned specimen set changed')
    barcodes = pd.read_csv(paths['barcodes'], header=None, sep='\t', dtype=str)[0]
    if not barcodes.is_unique or not assignments.index.isin(barcodes).all():
        raise ValueError('Hashing/matrix cell identities differ')
    obs = pd.DataFrame(index=barcodes)
    obs['source_tissue'] = 'unassigned'
    obs['source_sample'] = 'unassigned'; obs['source_donor'] = 'unassigned'
    obs.loc[assignments.index, 'source_tissue'] = 'PBMC'
    obs.loc[assignments.index, 'source_sample'] = assignments.source_sample_id
    obs.loc[assignments.index, 'source_donor'] = assignments.donor_id
    for column in ('singlet_hypothesis_probability', 'sensitivity_singlet_probability', 'Classification', 'assignment_method'):
        obs[column] = assignments[column].reindex(obs.index)
    paths['metadata'] = destination.with_suffix('.derived_metadata.csv')
    obs.to_csv(paths['metadata'])
    conversion = convert(paths, destination, cfg)
    return dict(kind='recomputed_hashsolo_RNA_counts', assets=manifests, hashing=report, conversion=conversion), list(paths.values())
