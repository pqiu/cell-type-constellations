"""Read a combined count triplet with explicitly reviewed specimen prefixes."""
import json
from pathlib import Path
import pandas as pd
from scripts.download_pbmc_asset import ROOT, download
from src.geo_mtx_table import convert


def prepare(cfg, destination):
    paths = {}; manifests = {}
    folder = ROOT / 'data/downloads'
    if cfg.get('asset_subdirectory'):
        name = cfg['asset_subdirectory']
        if Path(name).name != name or name in ['.', '..']:
            raise ValueError('Asset subdirectory must be a basename')
        folder = folder / name
    for role, asset in cfg['matrix_assets'].items():
        name = asset['filename']
        if role not in ('matrix', 'barcodes', 'features') or Path(name).name != name:
            raise ValueError('Unexpected count asset')
        paths[role] = download(asset['url'], folder / name, asset.get('expected_bytes'))
        manifests[role] = json.loads(paths[role].with_suffix(paths[role].suffix+'.manifest.json').read_text())
        if asset.get('sha256') and manifests[role]['sha256'] != asset['sha256']:
            raise ValueError('Reviewed original count asset checksum changed')
    barcodes = pd.read_csv(paths['barcodes'], sep='\t', header=None, dtype=str)[0]
    mode = cfg.get('barcode_key_mode', 'prefix')
    if mode == 'original_record_suffix':
        if cfg['study_id']!='GSE169246':raise ValueError('Original record-suffix selection requires source-specific review')
        if not barcodes.is_unique or not barcodes.str.contains('.',regex=False).all():
            raise ValueError('Invalid original record suffix barcodes')
        keys=barcodes.str.rsplit('.',n=1).str[-1]
        records=cfg['original_record_mapping']
        mapping=pd.DataFrame(records).set_index('source_sample_id')
        if not mapping.index.is_unique or set(keys)!=set(mapping.index):
            raise ValueError('Original RNA record set differs from exact barcode suffixes')
        for column in ['donor_id','source_tissue','source_GEO_accession']:
            if mapping[column].isna().any() or mapping[column].astype(str).str.strip().eq('').any():
                raise ValueError('Missing original record identity')
        obs=pd.DataFrame(dict(source_sample=keys.to_numpy(),
            source_donor=keys.map(mapping.donor_id).to_numpy(),
            source_tissue=keys.map(mapping.source_tissue).to_numpy(),
            source_GEO_accession=keys.map(mapping.source_GEO_accession).to_numpy()),index=barcodes)
    else:
        delimiter = '_' if mode == 'prefix' else '-'
        if mode not in ('prefix', 'suffix') or not barcodes.is_unique or not barcodes.str.contains(delimiter, regex=False).all():
            raise ValueError('Invalid source barcodes')
        keys = barcodes.str.split('_', n=1).str[0] if mode == 'prefix' else '-' + barcodes.str.rsplit('-', n=1).str[-1]
        expected = {str(row['source_sample_id']) for row in cfg['specimen_metadata']}
        if set(keys) != expected:
            raise ValueError('Barcode specimen prefixes differ from reviewed source metadata')
        donor_map = {str(row['source_sample_id']): str(row['donor_id']) for row in cfg['specimen_metadata']}
        if len(donor_map) != len(cfg['specimen_metadata']):
            raise ValueError('Duplicate specimen mapping')
        obs = pd.DataFrame({'source_sample': keys.to_numpy(), 'source_donor': keys.map(donor_map).to_numpy(),
                            'source_tissue': 'PBMC'}, index=barcodes)
    paths['metadata'] = destination.with_suffix('.derived_metadata.csv')
    obs.to_csv(paths['metadata'])
    report = convert(paths, destination, cfg)
    report['barcode_mapping'] = 'Exact reviewed ' + mode + '; complete source specimen set verified.'
    return {'kind': 'reviewed_specimen_prefix_triplet', 'assets': manifests, 'conversion': report}, list(paths.values())
