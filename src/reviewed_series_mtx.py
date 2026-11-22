"""Join a published cell/batch table to explicitly reviewed original specimens."""
import json
from pathlib import Path

import pandas as pd

from scripts.download_pbmc_asset import ROOT, download
from src.geo_mtx_table import convert


def original_assignments(path, cfg):
    original = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)
    if original.columns.tolist() != ['Cell', 'batch']:
        raise ValueError('Unexpected original cell-assignment schema')
    if not original.Cell.is_unique or original.eq('').any().any():
        raise ValueError('Ambiguous or missing original cell assignments')
    records = cfg['specimen_metadata']
    if len({r['source_sample_id'] for r in records}) != len(records):
        raise ValueError('Duplicate original batch mapping')
    if set(original.batch) != {r['source_sample_id'] for r in records}:
        raise ValueError('Original batch set differs from reviewed specimen mapping')
    donor = {r['source_sample_id']: r['donor_id'] for r in records}
    obs = pd.DataFrame(dict(source_original_cell_id=original.Cell.to_numpy(),
        source_original_batch=original.batch.to_numpy(), source_sample=original.batch.to_numpy(),
        source_donor=original.batch.map(donor).to_numpy(), source_tissue='PBMC'), index=original.Cell)
    if obs.source_donor.isna().any():
        raise ValueError('Unassigned donor in original batch table')
    return obs


def prepare(cfg, destination):
    from src.pbmc_pipeline import digest
    paths = {}; manifests = {}
    for role, asset in cfg['matrix_assets'].items():
        name = asset['filename']
        if role not in ['features', 'barcodes', 'matrix', 'cell_batch'] or Path(name).name != name:
            raise ValueError('Unexpected series count asset')
        p = download(asset['url'], ROOT / f'data/downloads/{cfg["dataset_id"]}_assets' / name, asset['expected_bytes'])
        if digest(p) != asset['sha256']:
            raise ValueError('Reviewed original series asset checksum changed')
        paths[role] = p
        manifests[role] = json.loads(p.with_suffix(p.suffix + '.manifest.json').read_text())
    obs = original_assignments(paths['cell_batch'], cfg)
    derived = destination.with_suffix('.derived_metadata.csv')
    obs.to_csv(derived)
    paths['metadata'] = derived
    conversion = convert(paths, destination, cfg)
    return dict(kind='original_series_counts_with_reviewed_batch_specimens', assets=manifests,
        conversion=conversion, original_sample_identity_review=cfg['source_metadata']['original_sample_identity_review']), list(paths.values())
