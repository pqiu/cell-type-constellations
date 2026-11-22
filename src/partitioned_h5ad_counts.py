"""Reunite disjoint cell partitions on exactly identical native RNA gene axes."""
import json
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem
from scipy import sparse

from scripts.download_pbmc_asset import ROOT, LIMIT, download, footprint


def prepare(cfg, destination):
    assets = cfg['source_assets']
    if len(assets) < 2 or not cfg.get('partition_reassembly_evidence'):
        raise ValueError('Multiple explicit source partitions and evidence required')
    matrices, observations, paths, manifests = [], [], [], []
    reference = None
    for asset in assets:
        name = asset['filename']
        if Path(name).name != name:
            raise ValueError('Unsafe source asset name')
        path = download(asset['url'], ROOT / 'data/downloads' / name, asset['expected_bytes'])
        manifest = json.loads(path.with_suffix(path.suffix + '.manifest.json').read_text())
        if asset.get('sha256') and manifest['sha256'] != asset['sha256']:
            raise ValueError('Pinned source partition checksum changed')
        if footprint() + 2 * sum(a['expected_bytes'] for a in assets) + 2_000_000_000 >= LIMIT:
            raise RuntimeError('Insufficient project disk budget for partition assembly')
        with h5py.File(path, 'r') as source:
            obs, var = read_elem(source['obs']), read_elem(source['raw/var'])
            if source['raw/X'].attrs['encoding-type'] != 'csr_matrix':
                raise ValueError('Explicit CSR native count source required')
            x = read_elem(source['raw/X'])
        if len(obs) != asset['expected_cells'] or x.shape != (len(obs), len(var)):
            raise ValueError('Source partition dimensions changed')
        if not obs.index.is_unique or not var.index.is_unique or (var.index.astype(str).str.strip() == '').any():
            raise ValueError('Invalid original cell/gene identifiers')
        if reference is None:
            reference = var.copy()
        elif not reference.index.equals(var.index) or not reference[cfg['source_gene_symbol_column']].equals(var[cfg['source_gene_symbol_column']]):
            raise ValueError('Partition gene axes differ; implicit alignment prohibited')
        x.check_format(full_check=True)
        if not x.has_canonical_format or not np.isfinite(x.data).all() or (x.data < 0).any() or not np.equal(x.data, np.floor(x.data)).all():
            raise ValueError('Invalid native integer RNA counts or repeated coordinates')
        if 'source_partition' in obs:
            raise ValueError('Refuse source metadata overwrite')
        obs['source_partition'] = asset['partition_label']
        observations.append(obs)
        matrices.append(x)
        paths.append(path)
        manifests.append(manifest)
    obs = pd.concat(observations, axis=0)
    if not obs.index.is_unique:
        raise ValueError('Source partitions overlap in cell identities; reassembly prohibited')
    counts = sparse.vstack(matrices, format='csr')
    staged = ad.AnnData(X=counts, obs=obs, var=reference)
    staged.write_h5ad(destination, compression='gzip')
    return dict(kind='disjoint_native_count_partition_reassembly', source_assets=manifests,
                partition_cells=[len(t) for t in observations], cells=len(obs), native_genes=len(reference),
                original_cell_ids_disjoint=True, identical_gene_axes=True, all_counts_finite_nonnegative_integer=True,
                arithmetic='Unchanged CSR rows concatenated, without expression transformation, feature selection or averaging.',
                partition_reassembly_evidence=cfg['partition_reassembly_evidence']), paths
