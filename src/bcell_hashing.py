"""Replay native paired RNA/CSP inference with pinned source and assignment IDs."""
import copy
import json
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scripts.download_pbmc_asset import ROOT
from scripts.audit_bcell_hashsolo_207475 import audit
from src.geo_assets import prepare_assets
from src.geo_inputs import convert_geo_tar


def prepare(cfg, destination):
    if cfg['study_id'] != 'GSE207475': raise ValueError('Reviewed source study required')
    archive, manifest, paths = prepare_assets(cfg)
    tables, report = audit()
    for key, value in cfg['hashing_settings'].items():
        if report[key] != value: raise ValueError('Pinned hashing settings changed: ' + key)
    pool_reports = {r['RNA_accession']:r for r in report['pools']}
    derived = copy.deepcopy(cfg); overrides = {}; specimens = set()
    for spec in derived['geo_samples']:
        gsm = spec['accession']; review = pool_reports[gsm]
        for field in ['assignment_sha256','source_count_sha256','source_tag_count_sha256']:
            if review[field] != spec[field]: raise ValueError('Pinned source/assignment checksum changed: ' + field)
        selected = tables[gsm]; specimens.update(selected.specimen_id)
        member = gsm + '_derived_hashsolo_assignments.csv'
        overrides[member] = selected.to_csv().encode()
        spec['cell_assignment'] = {'member':member,'barcode_prefix':'','donor_column':'biological_donor',
                                  'expected_donors':sorted(selected.biological_donor.unique())}
    if sorted(specimens) != cfg['expected_specimens']: raise ValueError('Physical participant set changed')
    convert_geo_tar(archive, destination, derived, metadata_overrides=overrides)
    converted = ad.read_h5ad(destination,backed='r')
    try:
        obs = converted.obs
        if set(obs['source_cell_specimen_id'].astype(str)) != specimens or set(obs.donor_id.astype(str)) != specimens:
            raise ValueError('Converted sample/donor identities differ')
        if not obs.groupby('source_cell_specimen_id',observed=True).donor_id.nunique().eq(1).all():
            raise ValueError('Several donors assigned to one specimen')
        if converted.n_vars != 36601 or converted.n_obs != cfg['expected_assigned_cells']:
            raise ValueError('Original RNA gene/cell coverage changed')
        for start in range(0,converted.n_obs,1024):
            stop=min(start+1024,converted.n_obs);x=sparse.csr_matrix(converted.X[start:stop]);x.eliminate_zeros()
            totals=pd.to_numeric(obs['source_cell_source_RNA_counts'].iloc[start:stop].astype(str)).to_numpy()
            genes=pd.to_numeric(obs['source_cell_source_RNA_genes'].iloc[start:stop].astype(str)).to_numpy()
            if not np.array_equal(np.asarray(x.sum(axis=1)).ravel(),totals) or not np.array_equal(np.diff(x.indptr),genes):
                raise ValueError('Converted RNA totals/support differ from original source audit')
    finally:
        converted.file.close()
    return {'kind':'native_RNA_and_observed_CSP_with_recomputed_hashsolo_assignments','assets':manifest,
            'hashing':report,'uncertainty':'Conditional sample-tag model assignments; unobserved tags are never zero-filled'},paths+[archive]
