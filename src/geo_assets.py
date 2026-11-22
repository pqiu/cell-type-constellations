"""Selectively retrieve count assets without downloading other assay modalities."""
import json
import shutil
import tarfile
from datetime import datetime, timezone
from pathlib import Path

from scripts.download_pbmc_asset import download, footprint, LIMIT, ROOT

def prepare_assets(cfg):
    folder = (ROOT / 'data' / 'downloads' / (cfg['study_id'] + '_assets')).resolve()
    download_root = (ROOT / 'data' / 'downloads').resolve()
    if not folder.is_relative_to(download_root): raise ValueError('Unsafe asset folder')
    paths = {}; manifests = []
    for sample in cfg['geo_samples']:
        for asset in sample['assets']:
            filename = asset['filename']
            if Path(filename).name != filename: raise ValueError('Asset name must be a basename')
            destination = download(asset['url'], folder / filename, asset.get('expected_bytes'))
            paths[filename] = destination
            manifests.append(json.loads(destination.with_suffix(destination.suffix + '.manifest.json').read_text()))
    archive = download_root / (cfg['dataset_id'] + '_SELECTED.tar')
    required = sum(p.stat().st_size for p in paths.values()) + 1024 * 1024
    if footprint() + required + 2_000_000_000 >= LIMIT or shutil.disk_usage(ROOT).free < required + 2_000_000_000:
        raise RuntimeError('Selected asset staging would exceed available disk budget')
    part = archive.with_suffix('.tar.part')
    with tarfile.open(part, 'w') as handle:
        for name, path in paths.items(): handle.add(path, arcname=name, recursive=False)
    part.replace(archive)
    manifest = dict(kind='selected_source_assets', retrieved_at=datetime.now(timezone.utc).isoformat(),
        study_id=cfg['study_id'], dataset_id=cfg['dataset_id'], assets=manifests,
        note='Local staging tar contains unchanged source bytes; it is not claimed to be an upstream archive.')
    return archive, manifest, list(paths.values())
