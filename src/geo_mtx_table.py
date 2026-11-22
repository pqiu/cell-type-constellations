"""Import a reviewed combined GEO matrix with barcode-indexed cell metadata."""
import gzip
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.io import mmread
from threadpoolctl import threadpool_limits

from scripts.download_pbmc_asset import ROOT, LIMIT, footprint, download


def convert(paths, destination, cfg):
    features = pd.read_csv(paths['features'], sep='\t', header=None, dtype=str, keep_default_na=False)
    barcodes = pd.read_csv(paths['barcodes'], sep='\t', header=None, dtype=str, keep_default_na=False)[0]
    if cfg.get('source_feature_format') == 'single_column_native_gene_symbols':
        if features.shape[1]!=1 or not features[0].is_unique or features[0].str.strip().eq('').any():
            raise ValueError('Unique nonempty original single-column gene symbols required')
        features[1]=features[0];features[2]='Gene Expression'
    obs = pd.read_csv(paths['metadata'], sep=cfg['metadata_separator'], index_col=cfg.get('metadata_index_column', 0))
    obs.index = obs.index.astype(str)
    if not obs.index.is_unique or not barcodes.is_unique:
        raise ValueError('Cell identities must be unique')
    if set(obs.index) != set(barcodes):
        raise ValueError('Metadata and matrix barcode sets differ')
    obs = obs.loc[barcodes].copy()
    if 'source_pbmc_values' in cfg:
        allowed = cfg['source_pbmc_values']
        if not isinstance(allowed, list) or not allowed or len(set(allowed)) != len(allowed):
            raise ValueError('Reviewed source selection values must be a nonempty unique list')
        selected = obs[cfg['source_tissue_column']].isin(allowed).to_numpy()
    else:
        selected = obs[cfg['source_tissue_column']].eq(cfg['source_pbmc_value']).to_numpy()
    if not selected.any():
        raise ValueError('No reviewed PBMC cells')
    permitted = set(cfg.get('allowed_feature_types', ['Gene Expression']))
    if features.shape[1] != 3 or not set(features[2]) <= permitted or 'Gene Expression' not in permitted:
        raise ValueError('Unexpected combined-matrix feature types')
    rna = features[2].eq('Gene Expression').to_numpy()
    if not rna.any():
        raise ValueError('No RNA features')
    if not features[0].is_unique:
        raise ValueError('Source feature identifiers are not unique')
    with gzip.open(paths['matrix'], 'rb') as handle, threadpool_limits(limits=8):
        x = mmread(handle, spmatrix=True)
    if x.shape != (len(features), len(barcodes)):
        raise ValueError('Matrix dimensions do not match genes and barcodes')
    if not np.isfinite(x.data).all() or (x.data<0).any():
        raise ValueError('Invalid original counts, including excluded tissue/features')
    if not np.issubdtype(x.data.dtype,np.integer) and not np.equal(x.data,np.floor(x.data)).all():
        raise ValueError('Original matrix must contain integer counts')
    x = x.T.tocsr()[selected][:, rna].copy()
    x.sum_duplicates();x.sort_indices();x.eliminate_zeros()
    features = features.loc[rna].copy()
    obs = obs.loc[selected].copy()
    mapping = pd.DataFrame(cfg['specimen_metadata']).set_index('source_sample_id')
    if not mapping.index.is_unique:
        raise ValueError('Specimen mapping is not unique')
    keys = obs[cfg['source_sample_column']].astype(str)
    if set(keys) != set(mapping.index):
        raise ValueError('Observed PBMC specimens differ from reviewed mapping')
    mapped = mapping.loc[keys].copy(); mapped.index = obs.index
    for column in mapped:
        if column in obs and not obs[column].astype(str).equals(mapped[column].astype(str)):
            raise ValueError('Source metadata conflict: ' + column)
        obs[column] = mapped[column]
    if not obs['donor_id'].astype(str).equals(obs[cfg['source_donor_column']].astype(str)):
        raise ValueError('Source donor and specimen mapping disagree')
    var = pd.DataFrame({'gene_symbols': features[1].to_numpy(), 'feature_type': features[2].to_numpy()}, index=features[0].to_numpy())
    if footprint() + x.data.nbytes + x.indices.nbytes + x.indptr.nbytes + 2_000_000_000 >= LIMIT:
        raise RuntimeError('Converted input would exceed project disk limit')
    ad.AnnData(X=x, obs=obs, var=var).write_h5ad(destination, compression='gzip')
    return {'source_cells': len(barcodes), 'selected_pbmc_cells': len(obs), 'genes': len(var),
            'source_sample_count': int(keys.nunique()), 'conversion': 'Original counts; barcode-aligned PBMC rows only; no expression transformation.'}


def prepare(cfg, destination):
    folder = ROOT / 'data' / 'downloads' / (cfg['dataset_id'] + '_assets')
    paths = {}; manifests = {}
    for role, asset in cfg['matrix_assets'].items():
        name = asset['filename']
        if Path(name).name != name:
            raise ValueError('Asset name must be a basename')
        paths[role] = download(asset['url'], folder / name, asset.get('expected_bytes'))
        manifests[role] = json.loads(paths[role].with_suffix(paths[role].suffix + '.manifest.json').read_text())
    conversion = convert(paths, destination, cfg)
    return {'kind': 'combined_geo_matrix', 'assets': manifests, 'conversion': conversion}, list(paths.values())
