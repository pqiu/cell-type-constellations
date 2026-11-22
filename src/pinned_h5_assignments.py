"""Join pinned author cell assignments and retain per-capture gene availability."""
import gzip
import hashlib
import json
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

from scripts.download_pbmc_asset import ROOT, LIMIT, download, footprint
from src.geo_inputs import align_pool_metadata


def prepare(cfg, destination):
    parts, paths, manifests, reviews = [], [], [], {}
    library_names = []
    symbol_by_id = {}
    for number, spec in enumerate(cfg['geo_samples']):
        asset, = spec['assets']
        name = asset['filename']
        if Path(name).name != name or name != spec['matrix_h5_member']:
            raise ValueError('Exact original H5 basename required')
        path = download(asset['url'], ROOT / f'data/downloads/{cfg["study_id"]}_assets' / name, asset.get('expected_bytes'))
        manifest = json.loads(path.with_suffix(path.suffix + '.manifest.json').read_text())
        if asset.get('sha256') and manifest['sha256'] != asset['sha256']:
            raise ValueError('Pinned original RNA source checksum changed')
        assignment = spec['cell_assignment']
        retained = assignment['retained_source']
        metadata_path = (ROOT / retained['path']).resolve()
        if not metadata_path.is_relative_to((ROOT / 'references/datasets').resolve()):
            raise ValueError('Assignments must be retained project metadata')
        raw = metadata_path.read_bytes()
        if hashlib.sha256(raw).hexdigest() != retained['sha256']:
            raise ValueError('Pinned author assignments changed')
        with h5py.File(path, 'r') as source:
            m, features = source['matrix'], source['matrix/features']
            ids, symbols = features['id'].asstr()[()], features['name'].asstr()[()]
            barcodes = pd.Series(m['barcodes'].asstr()[()])
            if set(features['feature_type'].asstr()[()]) != {'Gene Expression'}:
                raise ValueError('Only explicitly RNA-only source H5 files supported')
            values = m['data'][()]
            if not np.isfinite(values).all() or (values < 0).any() or not np.equal(values, np.floor(values)).all():
                raise ValueError('Invalid original integer RNA counts')
            counts = sparse.csc_matrix((values, m['indices'][()], m['indptr'][()]), shape=tuple(m['shape'][()])).T.tocsr()
        counts.check_format(full_check=True)
        if counts.shape != (len(barcodes), len(ids)) or not pd.Index(ids).is_unique or any(not v.strip() for v in ids):
            raise ValueError('Invalid native RNA gene/cell coordinates')
        counts.sum_duplicates()
        counts.sort_indices()
        counts.eliminate_zeros()
        gene_sha = hashlib.sha256(json.dumps([ids.tolist(), symbols.tolist()], separators=(',', ':')).encode()).hexdigest()
        prior = spec['source_count_review']
        if gene_sha != prior['gene_coordinate_sha256'] or list(counts.shape) != prior['shape']:
            raise ValueError('Original count dimensions or gene axis changed')
        for key, symbol in zip(ids, symbols):
            if key in symbol_by_id and symbol_by_id[key] != symbol:
                raise ValueError('Source libraries disagree on gene-ID-to-symbol mapping')
            symbol_by_id[key] = symbol
        decoded = gzip.decompress(raw) if raw.startswith(b'\x1f\x8b') else raw
        keep, assigned, review = align_pool_metadata(decoded, barcodes, assignment)
        counts = counts[keep]
        obs = pd.DataFrame(index=spec['accession'] + ':' + barcodes[keep].astype(str).to_numpy())
        obs.index.name = 'cell_id'
        obs['donor_id'] = assigned[assignment['donor_column']].to_numpy()
        obs['source_sample_id'] = obs.donor_id
        obs['source_library_id'] = spec['accession']
        obs['source_library_index'] = number
        obs['tissue'] = 'blood'
        obs['is_primary_data'] = True
        for column, value in spec['metadata'].items():
            if column in obs:
                raise ValueError('Refuse generated identity metadata overwrite')
            obs[column] = np.nan if value is None else value
        for column in assigned:
            obs['source_cell_' + column.replace('.', '_')] = pd.Categorical(assigned[column])
        parts.append((counts, pd.Index(ids), obs))
        paths.append(path)
        manifests.append(manifest)
        reviews[spec['accession']] = dict(**review, native_genes=len(ids), gene_coordinate_sha256=gene_sha,
                                          assignment_source=retained, source_sha256=manifest['sha256'])
        library_names.append(spec['accession'])
    if not parts or len(set(library_names)) != len(parts):
        raise ValueError('Nonempty unique source libraries required')
    union = pd.Index(sorted(symbol_by_id))
    availability = np.zeros((len(parts), len(union)), dtype=bool)
    matrices = []
    for number, (counts, ids, obs) in enumerate(parts):
        columns = union.get_indexer(ids)
        if (columns < 0).any():
            raise ValueError('Invalid native union mapping')
        availability[number, columns] = True
        coo = counts.tocoo()
        matrices.append(sparse.csr_matrix((coo.data, (coo.row, columns[coo.col])), shape=(len(obs), len(union))))
    combined = sparse.vstack(matrices, format='csr')
    obs = pd.concat([p[2] for p in parts])
    if not obs.index.is_unique:
        raise ValueError('Duplicate source-library cell identities')
    matrix_bytes = combined.data.nbytes + combined.indices.nbytes + combined.indptr.nbytes
    if footprint() + 3 * matrix_bytes + 2_000_000_000 >= LIMIT:
        raise RuntimeError('Assembled RNA input would exceed project disk budget')
    var = pd.DataFrame({'gene_symbols': [symbol_by_id[v] for v in union]}, index=union)
    output = ad.AnnData(X=combined, obs=obs, var=var)
    output.uns['source_gene_availability'] = availability
    output.uns['source_gene_ids'] = union.to_numpy()
    output.uns['source_library_names'] = library_names
    output.uns['source_assignment_audit'] = {library: {key: value for key, value in review.items()
        if isinstance(value, (int, np.integer))} for library, review in reviews.items()}
    output.write_h5ad(destination, compression='gzip')
    return dict(kind='pinned_author_assignments_with_native_H5_gene_availability', source_assets=manifests,
                source_libraries=reviews, native_gene_union=len(union),
                conversion='Native integer count values retained. Union zeros are storage placeholders with explicit per-library availability; incomplete all-cell gene statistics remain NaN, never observed-cell-only averages.'), paths
