"""Stage a verified native Matrix Market derivative with explicit capture metadata."""
import hashlib
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from scripts.download_pbmc_asset import ROOT, LIMIT, footprint, download


def pin_path(pin):
    path = (ROOT / pin['path']).resolve()
    if not path.is_relative_to(ROOT.resolve()) or not path.is_file():
        raise ValueError('Missing or unsafe pinned evidence')
    with path.open('rb') as stream:
        sha = hashlib.file_digest(stream, 'sha256').hexdigest()
    if sha != pin['sha256']:
        raise ValueError('Pinned evidence changed: ' + pin['path'])
    return path


def attach_capture_metadata(obs, mapping):
    if not obs.index.is_unique or not mapping.index.is_unique or set(obs.sampleID.astype(str)) != set(mapping.index):
        raise ValueError('Capture mapping must exactly cover unique source sample IDs')
    if set(obs).intersection(mapping.columns):
        raise ValueError('Capture mapping would overwrite source metadata')
    mapped = mapping.loc[obs.sampleID.astype(str)].copy()
    mapped.index = obs.index
    return obs.join(mapped, validate='one_to_one')


def verify_canonical_arrays(matrix, pins):
    if not matrix.has_canonical_format:
        raise ValueError('Canonical native CSR counts required')
    for name in ['data', 'indices', 'indptr']:
        values = getattr(matrix, name)
        if len(values) != pins[name]['elements']:
            raise ValueError('Native CSR array length changed')
        h = hashlib.sha256()
        for start in range(0, len(values), 3_000_000):
            block = values[start:start + 3_000_000]
            if name == 'data' and (not np.isfinite(block).all() or (block < 0).any() or not np.equal(block, np.floor(block)).all()):
                raise ValueError('Native RNA values must remain finite nonnegative integers')
            h.update(block.astype('<i8').tobytes())
        if h.hexdigest() != pins[name]['sha256']:
            raise ValueError('Native canonical CSR values changed')


def prepare(cfg, destination):
    audit_path = pin_path(cfg['native_count_audit'])
    audit = json.loads(audit_path.read_text(encoding='utf-8'))
    assert audit['complete'] and audit['all_original_stored_counts_finite_nonnegative_integer']
    original = audit['original_count_source']
    raw = download(original['url'], ROOT / original['path'], original['bytes'])
    pin_path(dict(path=raw.relative_to(ROOT).as_posix(), sha256=original['sha256']))
    bundle_path = ROOT / audit['native_candidate_bundle']['path']
    if not bundle_path.exists():
        from scripts.audit_PBMC_158055_original_counts import main as rebuild
        code = ROOT / 'scripts/audit_PBMC_158055_original_counts.py'
        if hashlib.sha256(code.read_bytes()).hexdigest() != audit['code_sha256']:
            raise ValueError('Pinned native candidate audit code changed')
        original_audit = audit_path.read_bytes()
        try:
            regenerated = rebuild()
            if {k: v for k, v in regenerated.items() if k != 'checked_at'} != {
                k: v for k, v in audit.items() if k != 'checked_at'}:
                raise ValueError('Rebuilt native candidate differs from immutable audit')
        finally:
            audit_path.write_bytes(original_audit)
    bundle = pin_path(audit['native_candidate_bundle'])
    source = ad.read_h5ad(bundle)
    assert source.shape == (340216, 27943) and source.obs_names.is_unique and source.var_names.is_unique
    verify_canonical_arrays(source.X, audit['candidate_raw_CSR_array_pins'])
    gene_table = pd.read_csv(pin_path(audit['native_gene_coordinates']), dtype=str, keep_default_na=False)
    assert source.var_names.tolist() == gene_table.native_gene_id.tolist()
    assert source.var.gene_symbols.astype(str).tolist() == gene_table.gene_symbols.tolist()
    identity = json.loads(pin_path(cfg['material_accounting_review']).read_text(encoding='utf-8'))
    assert identity['study_namespace'] == cfg['study_id']
    mapping = pd.read_csv(pin_path(identity['capture_metadata_map']), index_col='source_sample_id', dtype=str, keep_default_na=False)
    assert mapping.source_originality_as_reported.eq('Yes').all()
    source.obs = attach_capture_metadata(source.obs, mapping)
    source.obs['tissue'] = 'blood'
    source.obs['is_primary_data'] = True
    if footprint() + 2 * (source.X.data.nbytes + source.X.indices.nbytes + source.X.indptr.nbytes) + 2_000_000_000 >= LIMIT:
        raise RuntimeError('Native Matrix Market count staging exceeds project disk limit')
    provenance = dict(kind='pinned_native_MatrixMarket_candidate_bundle', native_count_audit=cfg['native_count_audit'],
                      original_count_source=original, native_candidate_bundle=audit['native_candidate_bundle'],
                      material_accounting_review=cfg['material_accounting_review'],
                      native_arrays_rechecked=True, source_gene_and_cell_axes_preserved=True,
                      conversion='Native original integer RNA counts and exact original capture metadata; no normalization, gene intersection or pooled capture centroids.')
    source.uns['source_preconversion_audit_json'] = json.dumps(provenance)
    source.write_h5ad(destination, compression='gzip')
    return provenance, [raw, bundle]
