"""Original HCA eQTL lane counts joined to deposited donor assignments."""
import hashlib
import io
import json
import zipfile

import anndata as ad
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.sparse import vstack
from threadpoolctl import threadpool_limits

from scripts.download_pbmc_asset import ROOT, LIMIT, download, footprint


def prepare(cfg, destination):
    folder = ROOT / 'data/downloads/hca_eqtl45'
    paths = {}; manifests = {}
    for role, asset in cfg['assets'].items():
        path = download(asset['url'], folder / asset['filename'], asset['expected_bytes'])
        with path.open('rb') as handle:
            if hashlib.file_digest(handle, 'sha256').hexdigest() != asset['sha256']:
                raise ValueError('HCA source checksum mismatch: ' + role)
        paths[role] = path
        manifests[role] = json.loads(path.with_suffix(path.suffix + '.manifest.json').read_text())
    with zipfile.ZipFile(paths['assignments']) as z:
        assignments = pd.read_csv(io.BytesIO(z.read('pilot3_persons.tsv')), sep='\t', header=None,
                                  names=['barcode', 'donor_id'], dtype=str).set_index('barcode')
    with zipfile.ZipFile(paths['labels']) as z:
        labels = pd.read_csv(io.BytesIO(z.read('barcodes_to_cell_types.tsv')), sep='\t', dtype=str).set_index('barcode')
    if not assignments.index.is_unique or not labels.index.is_unique or set(assignments.index) != set(labels.index):
        raise ValueError('Source assignment and label identities differ')
    donors = pd.DataFrame(cfg['donors']).set_index('donor_id')
    unknown = assignments.loc[~assignments.donor_id.isin(donors.index)]
    if unknown.donor_id.value_counts().to_dict() != cfg['excluded_assignment_counts']:
        raise ValueError('Unexpected unresolved donor assignments')
    matrices = []; observations = []; audits = {}; seen = set(); genes = None
    with zipfile.ZipFile(paths['counts']) as z, threadpool_limits(limits=8):
        for lane in range(1, 9):
            prefix = f'lane_{lane}/'
            features = pd.read_csv(io.BytesIO(z.read(prefix + 'genes.tsv')), sep='\t', header=None, dtype=str)
            if genes is None:
                genes = features
            if not features.equals(genes) or not features[0].is_unique:
                raise ValueError('Lane gene universes differ or repeat identifiers')
            barcodes = pd.read_csv(io.BytesIO(z.read(prefix + 'barcodes.tsv')), header=None, dtype=str)[0]
            ids = pd.Index([b.removesuffix('-1') + f'_lane{lane}' for b in barcodes])
            if not ids.is_unique or seen.intersection(ids):
                raise ValueError('Nonunique lane cell identity')
            seen.update(ids)
            x = mmread(io.BytesIO(z.read(prefix + 'matrix.mtx'))).T.tocsr()
            if x.shape != (len(ids), len(genes)) or not np.isfinite(x.data).all() or (x.data < 0).any() or (x.data != np.floor(x.data)).any():
                raise ValueError('Invalid original count matrix')
            mapped = assignments.reindex(ids)
            keep = mapped.donor_id.isin(donors.index).to_numpy()
            obs = mapped.loc[keep].copy()
            obs['source_sample_id'] = obs.donor_id
            obs['source_library_id'] = f'lane_{lane}'
            obs['cell_type'] = labels.reindex(obs.index).cell_type
            for column in donors:
                obs[column] = donors.loc[obs.donor_id, column].to_numpy()
            audits[f'lane_{lane}'] = dict(source_cells=len(ids), retained_assignment_cells=int(keep.sum()),
                                         unassigned_cells=int(mapped.donor_id.isna().sum()),
                                         unresolved_donor_cells=int((mapped.donor_id.notna() & ~mapped.donor_id.isin(donors.index)).sum()))
            matrices.append(x[keep]); observations.append(obs)
    if not set(assignments.index).issubset(seen):
        raise ValueError('Assigned cells missing from original matrices')
    obs = pd.concat(observations)
    if set(obs.donor_id) != set(donors.index):
        raise ValueError('Reviewed donors absent from counts')
    x = vstack(matrices, format='csr')
    if footprint() + x.data.nbytes + x.indices.nbytes + 2_000_000_000 >= LIMIT:
        raise RuntimeError('Conversion would exceed disk budget')
    var = pd.DataFrame({'gene_symbols': genes[1].to_numpy()}, index=genes[0].to_numpy())
    ad.AnnData(X=x, obs=obs, var=var, uns={'source_assignment_audit': audits}).write_h5ad(destination, compression='gzip')
    return {'kind': 'HCA original eQTL lane matrices', 'assets': manifests, 'assignment_audit': audits}, list(paths.values())
