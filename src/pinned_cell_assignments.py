"""Join retained, checksum-pinned source cell assignments to original counts."""
import hashlib
from pathlib import Path
from scripts.download_pbmc_asset import ROOT
from src.geo_assets import prepare_assets
from src.geo_inputs import convert_geo_tar


def prepare(cfg, destination):
    overrides = {}
    provenance = []
    for sample in cfg['geo_samples']:
        assignment = sample['cell_assignment']
        source = assignment['retained_source']
        path = (ROOT / source['path']).resolve()
        if not path.is_relative_to((ROOT / 'references/datasets').resolve()):
            raise ValueError('Assignment must be retained project metadata')
        raw = path.read_bytes()
        if hashlib.sha256(raw).hexdigest() != source['sha256']:
            raise ValueError('Source assignment checksum changed')
        member = assignment['member']
        if member in overrides:
            raise ValueError('Duplicate assignment member')
        overrides[member] = raw
        provenance.append(source)
    archive, manifest, paths = prepare_assets(cfg)
    convert_geo_tar(archive, destination, cfg, metadata_overrides=overrides)
    manifest['retained_cell_assignments'] = provenance
    return manifest, paths + [archive]
