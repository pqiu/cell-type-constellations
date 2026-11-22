"""Use pinned author cell assignments with original, unpadded GEO counts."""
import copy
import hashlib
import json
from pathlib import Path
import pandas as pd
from src.geo_inputs import convert_geo_tar
from scripts.download_pbmc_asset import ROOT, download


def assignments(cfg):
    asset = cfg['assignment_metadata']
    path = (ROOT / asset['path']).resolve()
    if not path.is_relative_to(ROOT.resolve()): raise ValueError('Metadata outside project')
    if hashlib.sha256(path.read_bytes()).hexdigest() != asset['sha256']:
        raise ValueError('Pinned author assignment metadata changed')
    frame = pd.read_csv(path, dtype=str)
    prefix = cfg.get('atlas_study_prefix', 'SCGT04_')
    if not frame.sampleID.str.startswith(prefix).all() or not frame.sampleID.str.endswith('_T0').all():
        raise ValueError('Expected reviewed baseline specimen identifiers')
    if set(frame.sampleID) != set(cfg['specimen_metadata']): raise ValueError('Clinical specimen set differs')
    if not frame.cell_id.is_unique: raise ValueError('Duplicate author cell IDs')
    used = set(); overrides = {}; derived = copy.deepcopy(cfg)
    for spec in derived['geo_samples']:
        rule = spec['atlas_selection']
        selected = frame.loc[frame.libraryID.eq(rule['library']) & frame.sampleID.isin(rule['samples'])].copy()
        if len(selected) != rule['expected_cells']: raise ValueError('Author selection changed')
        if used.intersection(selected.cell_id): raise ValueError('Cell assigned twice')
        used.update(selected.cell_id)
        barcode = selected.cell_id.str.rsplit('_', n=1).str[-1]
        expected_prefix = selected.libraryID + '_' + selected.sampleID.str.removeprefix(prefix) + '_'
        if not (expected_prefix + barcode).equals(selected.cell_id): raise ValueError('Cell identifier structure changed')
        if not barcode.str.fullmatch('[ACGT]{16}').all(): raise ValueError('Unexpected barcode')
        selected.index = barcode + '-1'
        if not selected.index.is_unique: raise ValueError('Duplicate barcode within library')
        selected['donor'] = selected.sampleID.str.removeprefix(prefix).str.removesuffix('_T0')
        for column in cfg.get('atlas_clinical_columns', ['age', 'sex', 'diseaseStatus']):
            selected[column] = selected.sampleID.map({s: m[column] for s,m in cfg['specimen_metadata'].items()})
        member = spec['accession'] + '_author_assignments.csv'
        overrides[member] = selected.to_csv().encode()
        spec['cell_assignment'] = dict(member=member, barcode_prefix='', donor_column='donor', expected_donors=sorted(selected.donor.unique()))
    if used != set(frame.cell_id): raise ValueError('Some retained author cells unassigned')
    return derived, overrides


def prepare(cfg, destination):
    derived, overrides = assignments(cfg)
    path = ROOT / 'data/downloads' / (cfg['dataset_id'] + '_RAW.tar')
    download(cfg['url'], path, cfg['expected_bytes'])
    manifest = json.loads(path.with_suffix('.tar.manifest.json').read_text())
    convert_geo_tar(path, destination, derived, metadata_overrides=overrides)
    return dict(kind='original_GEO_counts_with_pinned_atlas_assignments', counts=manifest,
                assignment_metadata=cfg['assignment_metadata']), [path]
