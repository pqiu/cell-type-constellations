"""Stage pinned native RNA cells with the authors' original pool/participant labels."""
import hashlib
import json

import anndata as ad
import numpy as np
import pandas as pd

from scripts.download_pbmc_asset import ROOT, LIMIT, footprint


def verify_arrays(matrix, pins):
    if not matrix.has_canonical_format: raise ValueError('Canonical native RNA arrays required')
    for name in ['data', 'indices', 'indptr']:
        values = getattr(matrix, name); assert len(values) == pins[name]['elements']
        h = hashlib.sha256()
        for start in range(0, len(values), 3_000_000):
            block = values[start:start + 3_000_000]
            if name == 'data': assert np.isfinite(block).all() and (block >= 0).all() and np.equal(block, np.floor(block)).all()
            h.update(block.astype('<i8').tobytes())
        assert h.hexdigest() == pins[name]['sha256']


def prepare(cfg, destination):
    from src.pbmc_pipeline import digest
    from scripts.prepare_PBMC_274284 import replay_sources
    audit = replay_sources()
    assert digest(ROOT / cfg['native_count_audit']['path']) == cfg['native_count_audit']['sha256']
    assert audit['complete'] and audit['exact_cross_capture_profile_pair_count'] == 0
    assert audit['all_original_feature_values_finite_nonnegative_integer']
    source = ROOT / audit['native_assigned_bundle']['path']; assert digest(source) == audit['native_assigned_bundle']['sha256']
    obj = ad.read_h5ad(source); verify_arrays(obj.X, audit['native_assigned_CSR_array_pins'])
    assert obj.shape == (audit['assigned_cell_observations'], 32738) and obj.obs_names.is_unique
    p = ROOT / cfg['source_identity_review']['path']; assert digest(p) == cfg['source_identity_review']['sha256']
    identity = json.loads(p.read_text(encoding='utf-8'))
    clinical = identity['source_clinical_records']
    assert set(obj.obs.donor_id.astype(str)).issubset(clinical)
    metadata = {
        'source_bulk_age_in_years_as_reported': 'age in_years', 'source_bulk_age_label_as_reported': 'age',
        'source_gender_as_reported': 'gender', 'source_condition_as_reported': 'condition',
        'source_phenotype_as_reported': 'phenotype', 'source_affected_gene_as_reported': 'gene'}
    for column, field in metadata.items():
        obj.obs[column] = pd.Categorical([clinical[code]['original_fields'].get(field, '') for code in obj.obs.donor_id.astype(str)])
    obj.obs['source_bulk_RNA_GEO_accession'] = pd.Categorical([clinical[code]['source_GSM'] for code in obj.obs.donor_id.astype(str)])
    obj.obs['source_capture_GEO_accession'] = pd.Categorical([identity['source_capture_records'][pool]['source_GSM'] for pool in obj.obs.source_pool.astype(str)])
    obj.obs['tissue'] = 'blood'; obj.obs['is_primary_data'] = True
    obj.obs['source_age_policy'] = 'Bulk-record ages are retained as reported; not assumed scRNA collection ages or converted clinical clocks.'
    obj.obs['source_assignment_policy'] = 'Published author SNG calls and BEST.LLK<-9.21034gate; LLK is not a calibrated posterior. No genotypes reconstructed.'
    obj.obs['source_replication_policy'] = 'Participant/material code groups;12dependent capture/enrichment views do not create independent donors. Exact visits/aliquots unspecified.'
    obj.obs['source_split_assignment'] = 'Unassigned; all captures/enrichment views of a source participant must share a split. Global aliases unresolved.'
    for column in obj.obs:
        if obj.obs[column].dtype == object: obj.obs[column] = pd.Categorical(obj.obs[column])
    needed = obj.X.data.nbytes + obj.X.indices.nbytes + obj.X.indptr.nbytes
    if footprint() + 2 * needed + 2_000_000_000 >= LIMIT: raise RuntimeError('Native RNA staging exceeds project disk budget')
    obj.write_h5ad(destination, compression='gzip')
    paths = [ROOT / audit['source_archive']['path'], source]
    assignment = json.loads((ROOT / audit['source_assignment_review']['path']).read_text(encoding='utf-8'))
    paths.append(ROOT / assignment['source_archive']['path'])
    paths.extend(ROOT / c['extracted_member']['path'] for c in audit['captures'])
    return dict(kind='original_RAG_PBMC_native_RNA', native_count_audit=cfg['native_count_audit'],
        source_identity_review=cfg['source_identity_review'], assigned_cells=obj.n_obs, native_genes=obj.n_vars,
        arithmetic='Unchanged author-assigned native RNA rows; no expression transformation or averaging. Antibody features excluded from RNA.',
        assignment_gate='Published author SNG and BEST.LLK<-9.21034'), paths
