"""Read a pinned native RNA export of the original GSE247917 source object."""
import hashlib
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import rdata

from scripts.download_pbmc_asset import ROOT, LIMIT, footprint, download
from src.rdata_counts import sparse_counts


def pinned_path(pin, root=ROOT):
    path = (root / pin['path']).resolve()
    if not path.is_relative_to(root.resolve()) or not path.is_file():
        raise ValueError('Missing or unsafe pinned source evidence')
    with path.open('rb') as stream:
        actual = hashlib.file_digest(stream, 'sha256').hexdigest()
    if actual != pin['sha256']:
        raise ValueError('Pinned source evidence changed: ' + pin['path'])
    return path


def verified_retention(plan, cells):
    if not plan.index.equals(cells) or not plan.index.is_unique:
        raise ValueError('Retention plan must exactly preserve ordered original cell IDs')
    if list(plan.columns) != ['retain_for_RNA_analysis'] or plan.iloc[:, 0].dtype != bool:
        raise ValueError('Retention decisions must be explicit nonmissing booleans')
    return plan.iloc[:, 0].to_numpy()


def verify_array_hash(values, dtype, pin):
    if len(values) != pin['elements'] or len(values) * np.dtype(dtype).itemsize != pin['bytes']:
        raise ValueError('Original sparse array length changed')
    h = hashlib.sha256()
    for start in range(0, len(values), 1_000_000):
        h.update(np.asarray(values[start:start + 1_000_000], dtype=dtype).tobytes())
    if h.hexdigest() != pin['sha256']:
        raise ValueError('Original sparse array values changed')


def prepare(cfg, destination):
    asset = cfg['source_asset']
    if Path(asset['filename']).name != asset['filename']:
        raise ValueError('Unsafe source filename')
    original = download(asset['url'], ROOT / 'data/downloads' / asset['filename'], asset['expected_bytes'])
    pinned_path(dict(path=original.relative_to(ROOT).as_posix(), sha256=asset['sha256']))
    export_path = pinned_path(cfg['native_export_receipt'])
    export = json.loads(export_path.read_text(encoding='utf-8'))
    assert export['source_object']['sha256'] == asset['sha256']
    validation_path = pinned_path(cfg['independent_native_validation'])
    validation = json.loads(validation_path.read_text(encoding='utf-8'))
    assert validation['original_export_receipt']['sha256'] == cfg['native_export_receipt']['sha256']
    repeat_path = pinned_path(cfg['exact_repeat_review'])
    repeat = json.loads(repeat_path.read_text(encoding='utf-8'))
    assert repeat['independent_native_validation']['sha256'] == cfg['independent_native_validation']['sha256']
    assert repeat['retention_plan'] == cfg['retention_plan']
    bundle = pinned_path(export['native_bundle'])
    source = rdata.read_rds(bundle)
    assert str(np.asarray(source['source_object_sha256']).item()) == asset['sha256']
    assert str(np.asarray(source['source_slot']).item()) == 'assays$RNA@counts'
    raw = source['counts']
    for name, array, dtype in [('data_float64', raw.x, '<f8'), ('indices_int32', raw.i, '<i4'),
                                ('indptr_int32', raw.p, '<i4')]:
        verify_array_hash(array, dtype, export['raw_CSC_array_pins'][name])
    counts, genes, cells = sparse_counts(raw)
    assert np.array_equal(genes.to_numpy(), source['native_gene_ids'])
    assert np.array_equal(cells.to_numpy(), source['native_cell_ids'])
    metadata = source['metadata']
    assert metadata.index.astype(str).equals(cells)
    original_metadata = rdata.read_rds(pinned_path(export['original_non_receptor_metadata_RDS']))
    pd.testing.assert_frame_equal(metadata, original_metadata, check_exact=True)
    plan = pd.read_csv(pinned_path(cfg['retention_plan']), index_col=0)
    keep = verified_retention(plan, cells)
    assert int(keep.sum()) == repeat['proposed_analysis_rows'] == 189771
    obs = metadata.loc[cells[keep]].copy()
    counts = counts[keep].copy()
    assignments = pd.read_csv(pinned_path(cfg['metadata_assignments']), index_col=0, dtype=str,
                              keep_default_na=False)
    if not assignments.index.equals(obs.index) or set(assignments).intersection(obs):
        raise ValueError('Source assignments must match exact cell order and add only new fields')
    for col in assignments:
        obs[col] = assignments[col]
    for col in obs:
        if obs[col].dtype == object or isinstance(obs[col].dtype, pd.StringDtype):
            obs[col] = obs[col].astype('string').fillna('not_reported').astype(str)
    obs['author_cell_type'] = obs['celltype.l2'].astype(str)
    obs['tissue'] = 'blood'
    obs['is_primary_data'] = True
    if footprint() + 2 * (counts.data.nbytes + counts.indices.nbytes + counts.indptr.nbytes) + 2_000_000_000 >= LIMIT:
        raise RuntimeError('Original native RNA staging exceeds project disk limit')
    var = pd.DataFrame({'gene_symbols': genes.to_numpy()}, index=genes)
    audit = dict(kind='verified_original_GSE247917_native_RNA', source_asset=asset,
                 native_export_receipt=cfg['native_export_receipt'],
                 independent_native_validation=cfg['independent_native_validation'],
                 exact_repeat_review=cfg['exact_repeat_review'], metadata_assignments=cfg['metadata_assignments'],
                 original_cells=len(cells), retained_original_observations=len(obs),
                 repeated_RNA_rows_excluded=int((~keep).sum()),
                 raw_sparse_array_hashes_rechecked=True, native_gene_axis_preserved=True,
                 conversion='Original RNA counts only, exact proved repeat map; no normalization, alignment, or imputation.',
                 reader_version=rdata.__version__)
    result = ad.AnnData(X=counts, obs=obs, var=var)
    result.uns['source_preconversion_audit_json'] = json.dumps(audit)
    result.write_h5ad(destination, compression='gzip')
    return audit, [original, bundle]
