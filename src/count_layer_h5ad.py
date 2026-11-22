"""Stage one explicitly audited count layer as backed X without loading it."""
import hashlib
import json
import shutil
from pathlib import Path

import h5py

from scripts.download_pbmc_asset import ROOT, LIMIT, download, footprint


def prepare(cfg, destination):
    asset = cfg['source_asset']
    name = asset['filename']
    if Path(name).name != name:
        raise ValueError('Source asset must be a basename')
    source = download(asset['url'],ROOT / 'data/downloads' / name,asset['expected_bytes'])
    audit_path = (ROOT / cfg['source_count_audit']['path']).resolve()
    if not audit_path.is_relative_to((ROOT / 'references/datasets').resolve()):
        raise ValueError('Source count audit must be retained project metadata')
    raw_audit = audit_path.read_bytes()
    if hashlib.sha256(raw_audit).hexdigest() != cfg['source_count_audit']['sha256']:
        raise ValueError('Source count audit changed')
    audit = json.loads(raw_audit)
    if not audit['complete'] or not audit['all_values_finite_nonnegative_integer'] or not audit['all_native_CSR_coordinates_checked'] or audit['duplicate_native_gene_entries']:
        raise ValueError('Full native count/coordinate audit required')
    if source.stat().st_size != audit['source_manifest']['bytes']:
        raise ValueError('Audited source file size changed')
    source_manifest = json.loads(source.with_suffix(source.suffix + '.manifest.json').read_text())
    if source_manifest['sha256'] != audit['source_manifest']['sha256']:
        raise ValueError('Audited source checksum changed')
    allowance = source.stat().st_size + 2_000_000_000
    if footprint() + allowance >= LIMIT or shutil.disk_usage(ROOT).free < allowance:
        raise RuntimeError('Insufficient disk budget for count-layer staging')
    layer = cfg['count_layer']
    if not layer or '/' in layer:
        raise ValueError('An explicit single layer name is required')
    partial = destination.with_suffix(destination.suffix + '.part')
    with h5py.File(source,'r') as original,h5py.File(partial,'w') as output:
        matrix = original['layers'][layer]
        if matrix.attrs['encoding-type'] != 'csr_matrix' or list(matrix.attrs['shape']) != audit['shape']:
            raise ValueError('Audited native count layer dimensions changed')
        output.attrs['encoding-type'] = 'anndata'
        output.attrs['encoding-version'] = '0.1.0'
        # HDF5 copies the stored groups without materializing a whole sparse
        # matrix. No source X, extra layers, graphs, or embeddings are staged.
        original.copy(matrix,output,name='X')
        original.copy('obs',output)
        original.copy('var',output)
    partial.replace(destination)
    return dict(kind='audited_native_count_layer',source=source_manifest,count_layer=layer,
        source_count_audit=cfg['source_count_audit'],staging='Unchanged native CSR group copied to backed X; original obs/var preserved. No count transformation.'),[source]
