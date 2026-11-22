"""Read reviewed sparse R count objects without executing serialized R code."""
import json
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import rdata
from scipy import sparse
from scripts.download_pbmc_asset import ROOT, LIMIT, footprint, download


def sparse_counts(obj):
    if list(getattr(obj, 'class')) != ['dgCMatrix']:
        raise ValueError('Expected reviewed dgCMatrix count object')
    shape = tuple(int(v) for v in obj.Dim)
    if len(shape) != 2 or len(obj.p) != shape[1] + 1 or obj.p[0] != 0 or obj.p[-1] != len(obj.x):
        raise ValueError('Invalid sparse column pointers or dimensions')
    if len(obj.i) != len(obj.x) or (np.diff(obj.p) < 0).any() or (obj.i < 0).any() or (obj.i >= shape[0]).any():
        raise ValueError('Invalid sparse indices')
    for start in range(0, len(obj.x), 1_000_000):
        values = obj.x[start:start+1_000_000]
        if not np.isfinite(values).all() or (values < 0).any() or not np.equal(values, np.floor(values)).all():
            raise ValueError('Source does not contain finite nonnegative integer counts')
    genes, cells = map(pd.Index, obj.Dimnames)
    if len(genes) != shape[0] or len(cells) != shape[1] or not genes.is_unique or not cells.is_unique:
        raise ValueError('Gene/cell identities must be unique and match dimensions')
    matrix = sparse.csc_matrix((obj.x, obj.i, obj.p), shape=shape)
    if not matrix.has_canonical_format:
        raise ValueError('Sparse matrix has duplicated or unsorted coordinates')
    return matrix.T.tocsr(), genes.astype(str), cells.astype(str)


def prepare(cfg, destination):
    paths = {}; manifests = {}
    for role, asset in cfg['matrix_assets'].items():
        name = asset['filename']
        if Path(name).name != name: raise ValueError('Unsafe source asset name')
        paths[role] = download(asset['url'], ROOT / 'data/downloads' / name, asset.get('expected_bytes'))
        manifests[role] = json.loads(paths[role].with_suffix(paths[role].suffix+'.manifest.json').read_text())
    objects = rdata.read_rda(paths['counts'])
    if set(objects) != {cfg['r_object_name']}: raise ValueError('Unexpected R object names')
    counts, genes, cells = sparse_counts(objects[cfg['r_object_name']])
    obs = pd.read_csv(paths['metadata'], index_col=0, dtype=str)
    if not obs.index.is_unique or set(obs.index) != set(cells): raise ValueError('Cell metadata does not exactly match count barcodes')
    obs = obs.loc[cells].copy()
    mapping = pd.DataFrame(cfg['specimen_metadata']).set_index('source_sample_id')
    keys = obs[cfg['source_sample_column']]
    if not mapping.index.is_unique or set(keys) != set(mapping.index): raise ValueError('Unreviewed sample identities')
    mapped = mapping.loc[keys].copy(); mapped.index = obs.index
    for col in mapped:
        if col in obs: raise ValueError('Mapping would overwrite source metadata')
        obs[col] = mapped[col]
    obs['source_sample_id'] = keys.to_numpy()
    obs['tissue'] = 'blood'; obs['is_primary_data'] = True
    if footprint() + 2*(counts.data.nbytes+counts.indices.nbytes+counts.indptr.nbytes) + 2_000_000_000 >= LIMIT:
        raise RuntimeError('Count conversion exceeds project disk limit')
    var = pd.DataFrame({'gene_symbols': genes.to_numpy()}, index=genes)
    ad.AnnData(X=counts, obs=obs, var=var).write_h5ad(destination, compression='gzip')
    return dict(kind='reviewed_rdata_counts', assets=manifests, reader_version=rdata.__version__,
        conversion='Exact dgCMatrix counts transposed to cells by genes; barcode metadata join; no normalization.'), list(paths.values())
