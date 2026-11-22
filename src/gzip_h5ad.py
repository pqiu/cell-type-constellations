"""Import a gzipped source H5AD with a reviewed library-to-specimen crosswalk."""
import gzip
import json
import shutil
from pathlib import Path
import anndata as ad
import pandas as pd
from scripts.download_pbmc_asset import ROOT, LIMIT, footprint, download


def prepare(cfg, destination):
    asset = cfg['source_asset']; name = asset['filename']
    if Path(name).name != name or not name.endswith('.h5ad.gz'): raise ValueError('Unsafe compressed H5AD source name')
    folder = ROOT / 'data/downloads'
    if cfg.get('source_asset_subdirectory'):
        subdirectory = cfg['source_asset_subdirectory']
        if Path(subdirectory).name != subdirectory or subdirectory in {'.', '..'}:
            raise ValueError('Unsafe compressed source subdirectory')
        folder = folder / subdirectory
    path = download(asset['url'], folder / name, asset['expected_bytes'])
    manifest = json.loads(path.with_suffix(path.suffix+'.manifest.json').read_text())
    if asset.get('sha256') and manifest.get('sha256') != asset['sha256']:
        raise ValueError('Pinned original source checksum changed')
    expanded = path.with_suffix('')
    allowance = min(LIMIT-footprint(), shutil.disk_usage(ROOT).free)-2_000_000_000
    total = 0
    with gzip.open(path,'rb') as source, expanded.open('wb') as target:
        while chunk := source.read(8*1024*1024):
            total += len(chunk)
            if total > allowance: raise RuntimeError('Gzip expansion exceeds disk budget')
            target.write(chunk)
    source = ad.read_h5ad(expanded)
    mapping = pd.DataFrame(cfg['specimen_metadata']).set_index('library')
    keys = source.obs[cfg['source_library_column']].astype(str)
    if not mapping.index.is_unique or set(keys) != set(mapping.index): raise ValueError('Unreviewed library identity')
    mapped = mapping.loc[keys].copy(); mapped.index = source.obs_names
    for column in mapped:
        if column in source.obs: raise ValueError('Refuse to overwrite source metadata')
        source.obs[column] = mapped[column]
    genes = (pd.Series(source.var_names.astype(str), index=source.var_names)
             if cfg['source_gene_symbol_column'] == '__index__'
             else source.var[cfg['source_gene_symbol_column']].astype(str))
    if not genes.is_unique or genes.str.strip().eq('').any(): raise ValueError('Invalid source gene coordinates')
    identifiers = genes
    if 'source_gene_id_column' in cfg:
        identifiers = source.var[cfg['source_gene_id_column']].astype(str)
        if not identifiers.is_unique or identifiers.str.strip().eq('').any():
            raise ValueError('Invalid original gene IDs')
    if 'expected_feature_types' in cfg:
        if set(source.var[cfg['feature_type_column']].astype(str)) != set(cfg['expected_feature_types']):
            raise ValueError('Source feature types changed')
    source.var['original_feature_index'] = source.var_names.astype(str)
    source.var_names = identifiers.to_numpy()
    source.var['gene_symbols'] = genes.to_numpy()
    for key in list(source.obsm): del source.obsm[key]
    for key in list(source.obsp): del source.obsp[key]
    source.uns.clear()
    if footprint()+total+2_000_000_000 >= LIMIT: raise RuntimeError('Converted H5AD exceeds disk budget')
    source.write_h5ad(destination,compression='gzip')
    coordinate_policy = ('Original gene IDs become coordinates; original symbols retained.'
        if 'source_gene_id_column' in cfg else 'Source gene symbols become coordinates.')
    return dict(kind='reviewed_gzip_h5ad',source=manifest,expanded_bytes=total,
        crosswalk=cfg['crosswalk_evidence'],conversion='Original X counts; exact library crosswalk; ' + coordinate_policy), [path,expanded]
