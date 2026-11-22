"""Preserve original Parse counts while joining sublibraries by clinical specimen."""
import gzip
import hashlib
import json
import shutil
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
from threadpoolctl import threadpool_limits

from scripts.download_pbmc_asset import ROOT, LIMIT, footprint, download


def sha(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def convert(gene_path, libraries, destination, cfg):
    """Write one library at a time; no dense array or full-matrix concatenation."""
    genes = pd.read_csv(gene_path, keep_default_na=False)
    if list(genes.columns) != ['gene_id', 'gene_name', 'genome'] or not genes.gene_id.is_unique:
        raise ValueError('Unreviewed Parse gene coordinates')
    if genes.gene_id.str.strip().eq('').any() or set(genes.genome) != {'hg38'}:
        raise ValueError('Missing original gene identity or unexpected reference')
    if sha(gene_path) != cfg['gene_asset']['sha256']:
        raise ValueError('Original gene axis changed')
    metadata = pd.DataFrame(cfg['specimen_metadata']).set_index('source_sample_id')
    if not metadata.index.is_unique or metadata.donor_id.isna().any():
        raise ValueError('Specimen and recipient crosswalk must be explicit')
    frames = []
    selections = []
    for spec, paths in libraries:
        if sha(paths['metadata']) != spec['metadata_asset']['sha256']:
            raise ValueError('Original cell metadata changed')
        obs = pd.read_csv(paths['metadata'], keep_default_na=False)
        if not obs.bc_wells.is_unique or obs.bc_wells.str.strip().eq('').any() or set(obs.species) != {'hg38'}:
            raise ValueError('Invalid Parse source cell coordinates')
        labels = obs['sample']
        allowed = set(metadata.index) | set(cfg['excluded_graft_samples'])
        if not set(labels).issubset(allowed):
            raise ValueError('Unexpected clinical sample assignment')
        keep = labels.isin(metadata.index).to_numpy()
        selected = obs.loc[keep].copy()
        selected.index = pd.Index(spec['accession'] + ':' + selected.bc_wells.astype(str))
        selected.index.name = None
        selected['source_library_id'] = spec['accession']
        selected['source_sample_id'] = selected['sample'].to_numpy()
        selected['tissue'] = 'blood'
        for column in metadata.columns:
            if column in selected:
                raise ValueError('Refuse to overwrite deposited cell metadata')
            selected[column] = metadata.loc[selected['source_sample_id'], column].to_numpy()
        frames.append(selected)
        selections.append((obs, keep))
    combined_obs = pd.concat(frames)
    if not combined_obs.index.is_unique or set(combined_obs.source_sample_id) != set(metadata.index):
        raise ValueError('Missing biological specimen or duplicate source cell identity')
    var = pd.DataFrame({'gene_symbols': genes.gene_name.to_numpy(),
                        'source_genome': genes.genome.to_numpy()}, index=genes.gene_id.to_numpy())
    part = Path(str(destination) + '.part')
    if part.exists():
        raise ValueError('Inspect prior derived conversion partial before retrying')
    if footprint() + 2_000_000_000 >= LIMIT:
        raise RuntimeError('Project disk limit would be exceeded')
    ad.AnnData(sparse.csr_matrix((len(combined_obs), len(var)), dtype=np.int64),
               obs=combined_obs, var=var).write_h5ad(part, compression='gzip')
    audits = {}
    row_offset = nnz_offset = total_umis = 0
    with h5py.File(part, 'r+') as handle:
        matrix = handle['X']
        for name in ['data', 'indices', 'indptr']:
            del matrix[name]
        data = matrix.create_dataset('data', shape=(0,), maxshape=(None,), dtype='<i8', compression='gzip')
        indices = matrix.create_dataset('indices', shape=(0,), maxshape=(None,), dtype='<i4', compression='gzip')
        indptr = matrix.create_dataset('indptr', shape=(len(combined_obs) + 1,), dtype='<i8', compression='gzip')
        indptr[0] = 0
        for (spec, paths), (obs, keep) in zip(libraries, selections):
            if sha(paths['matrix']) != spec['matrix_asset']['sha256']:
                raise ValueError('Original matrix differs from ordered-coordinate audit')
            with gzip.open(paths['matrix'], 'rb') as stream, threadpool_limits(8):
                x = mmread(stream, spmatrix=True).tocsr()
            x.check_format(full_check=True)
            if x.shape != (len(obs), len(var)) or not np.isfinite(x.data).all() or (x.data < 0).any() or not np.equal(x.data, np.floor(x.data)).all():
                raise ValueError('Invalid original Parse count values/dimensions')
            x.sum_duplicates(); x.sort_indices(); x.eliminate_zeros()
            if not np.array_equal(np.asarray(x.sum(axis=1)).ravel(), obs.tscp_count.to_numpy()) or not np.array_equal(np.diff(x.indptr), obs.gene_count.to_numpy()):
                raise ValueError('Original count rows disagree with deposited cell metadata')
            x = x[keep].tocsr()
            if x.data.size and x.data.max() > np.iinfo(np.int64).max:
                raise ValueError('Count exceeds exact signed integer storage')
            required = x.nnz * 12 + x.indptr.nbytes
            if footprint() + required + 2_000_000_000 >= LIMIT or shutil.disk_usage(ROOT).free < required + 2_000_000_000:
                raise RuntimeError('Streaming conversion would exceed available disk budget')
            end = nnz_offset + x.nnz
            data.resize((end,)); indices.resize((end,))
            data[nnz_offset:end] = x.data
            indices[nnz_offset:end] = x.indices
            indptr[row_offset + 1:row_offset + x.shape[0] + 1] = x.indptr[1:].astype(np.int64) + nnz_offset
            audits[spec['accession']] = {'source_cells': len(obs), 'selected_PBMC_cells': x.shape[0],
                                         'excluded_graft_cells': int((~keep).sum()), 'nonzeros': int(x.nnz)}
            row_offset += x.shape[0]; nnz_offset = end
            total_umis += int(x.sum())
            print({'converted_library': spec['accession'], 'PBMC_cells': x.shape[0]}, flush=True)
        if row_offset != len(combined_obs) or int(indptr[-1]) != nnz_offset:
            raise ValueError('Incomplete streamed CSR matrix')
    check = ad.read_h5ad(part, backed='r')
    try:
        if check.shape != (len(combined_obs), len(var)) or not check.obs.index.equals(combined_obs.index) or not check.var.index.equals(var.index):
            raise ValueError('Streamed AnnData coordinate round-trip failed')
    finally:
        check.file.close()
    part.replace(destination)
    return {'libraries': audits, 'PBMC_cells': len(combined_obs), 'source_features': len(var),
            'PBMC_UMIs': total_umis, 'biological_samples': combined_obs.source_sample_id.nunique(),
            'clinical_recipients': combined_obs.donor_id.nunique(),
            'conversion': 'Unchanged integer RNA counts; graft rows excluded, all original gene columns retained. Sublibraries join by explicit source visit and clinical recipient IDs.'}


def prepare(cfg, destination):
    folder = ROOT / 'data/downloads' / (cfg['study_id'] + '_assets')
    gene = cfg['gene_asset']
    gene_path = download(gene['url'], folder / gene['filename'])
    libraries = []
    manifests = {}
    counts = []
    for spec in cfg['parse_libraries']:
        paths = {}
        for role in ['metadata', 'matrix']:
            asset = spec[role + '_asset']
            if Path(asset['filename']).name != asset['filename']:
                raise ValueError('Asset filename must be a basename')
            path = download(asset['url'], folder / asset['filename'])
            paths[role] = path
            manifests[asset['filename']] = json.loads(path.with_suffix(path.suffix + '.manifest.json').read_text())
            if role == 'matrix':
                counts.append(path)
        libraries.append((spec, paths))
    conversion = convert(gene_path, libraries, destination, cfg)
    manifests[gene['filename']] = json.loads(gene_path.with_suffix(gene_path.suffix + '.manifest.json').read_text())
    return {'kind': 'original_Parse_GEO_RNA', 'assets': manifests, 'conversion': conversion}, counts
