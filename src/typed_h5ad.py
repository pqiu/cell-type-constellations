"""Select explicitly typed RNA features and reviewed cell identities from H5AD."""
import json
from pathlib import Path
import anndata as ad
import numpy as np
from scripts.download_pbmc_asset import ROOT, LIMIT, footprint, download
from src.source_selection import source_exclusion_mask


def prepare(cfg, destination):
    asset = cfg['source_asset']; name = asset['filename']
    if Path(name).name != name: raise ValueError('Unsafe source asset name')
    path = download(asset['url'], ROOT / 'data/downloads' / name, asset['expected_bytes'])
    manifest = json.loads(path.with_suffix(path.suffix+'.manifest.json').read_text())
    source = ad.read_h5ad(path, backed='r')
    try:
        features = source.var[cfg['feature_type_column']].eq('Gene Expression').to_numpy()
        excluded = source_exclusion_mask(source.obs, cfg['source_cell_exclusions'])
        if not features.any() or excluded.all(): raise ValueError('No eligible RNA observations')
        if footprint() + 2*path.stat().st_size + 2_000_000_000 >= LIMIT:
            raise RuntimeError('H5AD conversion exceeds project disk budget')
        selected = source[~excluded, features].to_memory()
        ids = selected.var[cfg['gene_id_column']].astype(str)
        if not ids.is_unique or ids.str.strip().eq('').any(): raise ValueError('Invalid source gene identifiers')
        selected.var['source_gene_symbols'] = selected.var_names.astype(str)
        selected.var['gene_symbols'] = selected.var_names.astype(str)
        selected.var_names = ids.to_numpy()
        for key in list(selected.obsm): del selected.obsm[key]
        for key in list(selected.obsp): del selected.obsp[key]
        selected.uns.clear()
        selected.write_h5ad(destination, compression='gzip')
        review = dict(source_cells=source.n_obs, retained_source_cells=selected.n_obs,
            excluded_source_cells=int(excluded.sum()), source_features=source.n_vars,
            retained_rna_features=int(features.sum()), excluded_feature_types=source.var.loc[~features,cfg['feature_type_column']].value_counts().to_dict(),
            gene_mapping='Source gene IDs become coordinates; original symbols retained. No count transformation.')
    finally:
        source.file.close()
    return dict(kind='typed_source_h5ad',source=manifest,conversion=review), [path]
