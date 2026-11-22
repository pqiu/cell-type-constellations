"""Replay pinned sample-tag inference, keeping RNA values and visit identities."""
import copy
import hashlib
import json
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scripts.download_pbmc_asset import ROOT
from scripts.audit_cll_hashsolo_283397 import audit
from src.geo_assets import prepare_assets
from src.geo_inputs import convert_geo_tar


def prepare(cfg, destination):
    if cfg['study_id'] != 'GSE283397': raise ValueError('Reviewed source study required')
    archive, manifest, paths = prepare_assets(cfg)
    tables, report = audit()
    for key in ['scanpy_version', 'hashsolo_code_sha256', 'priors', 'seed', 'number_of_noise_barcodes', 'tag_column_order', 'minimum_singlet_probability', 'RNA_QC']:
        if report[key] != cfg['hashing_settings'][key]: raise ValueError('Pinned hashing configuration changed: ' + key)
    pool_reports = {r['pool']: r for r in report['pools']}
    derived = copy.deepcopy(cfg); overrides = {}; specimen_ids = set(); donor_ids = set()
    for spec in derived['geo_samples']:
        pool = spec['source_pool']; review = pool_reports[pool]
        if review['assignment_sha256'] != spec['assignment_sha256'] or review['source_count_sha256'] != spec['source_count_sha256']:
            raise ValueError('Pinned source counts or recomputed assignments changed')
        selected = tables[pool]
        member = spec['accession'] + '_derived_hashsolo_assignments.csv'
        overrides[member] = selected.to_csv().encode()
        spec['cell_assignment'] = {'member': member, 'barcode_prefix': '', 'donor_column': 'biological_donor',
                                   'expected_donors': sorted(selected.biological_donor.unique())}
        specimen_ids.update(selected.specimen_id); donor_ids.update(selected.biological_donor)
    if sorted(specimen_ids) != cfg['expected_specimens'] or sorted(donor_ids) != cfg['expected_donors']:
        raise ValueError('Physical donor or specimen set changed')
    convert_geo_tar(archive, destination, derived, metadata_overrides=overrides)
    # Check the actual conversion against independently computed source totals
    # and support for every selected cell, including the repeated-visit joins.
    converted = ad.read_h5ad(destination, backed='r')
    try:
        obs = converted.obs
        if set(obs['source_cell_specimen_id'].astype(str)) != specimen_ids or set(obs.donor_id.astype(str)) != donor_ids:
            raise ValueError('Converted specimen/donor identities differ')
        if not obs.groupby('source_cell_specimen_id', observed=True).donor_id.nunique().eq(1).all():
            raise ValueError('Several donors assigned to one physical specimen')
        for start in range(0, converted.n_obs, 1024):
            stop = min(start + 1024, converted.n_obs)
            x = sparse.csr_matrix(converted.X[start:stop]); x.eliminate_zeros()
            expected_total = pd.to_numeric(obs['source_cell_source_RNA_counts'].iloc[start:stop].astype(str)).to_numpy()
            expected_genes = pd.to_numeric(obs['source_cell_source_RNA_genes'].iloc[start:stop].astype(str)).to_numpy()
            if not np.array_equal(np.asarray(x.sum(axis=1)).ravel(), expected_total) or not np.array_equal(np.diff(x.indptr), expected_genes):
                raise ValueError('Converted RNA counts/support differ from original matrix audit')
    finally:
        converted.file.close()
    return {'kind': 'original_RNA_with_recomputed_hashsolo_assignments', 'assets': manifest,
            'hashing': report, 'sample_key': 'source_cell_specimen_id; donor identity remains biological_donor'}, paths + [archive]
