"""Stage native RNA rows whose original authors explicitly identify each cell."""
import json

import anndata as ad
import pandas as pd

from scripts.download_pbmc_asset import ROOT, LIMIT, footprint
from src.rag_pbmc_original import verify_arrays


def prepare(cfg, destination):
    from src.pbmc_pipeline import digest
    from scripts.prepare_PBMC_266519 import replay_sources
    audit = replay_sources()
    assert digest(ROOT / cfg['native_count_audit']['path']) == cfg['native_count_audit']['sha256']
    assert audit['complete'] and audit['all_original_stored_values_finite_nonnegative_integer']
    bundle = ROOT / audit['native_selected_bundle']['path']; assert digest(bundle) == audit['native_selected_bundle']['sha256']
    obj = ad.read_h5ad(bundle); verify_arrays(obj.X, audit['native_selected_CSR_array_pins'])
    assert obj.shape == (89753, 36601) and obj.obs_names.is_unique
    ip = ROOT / cfg['source_identity_review']['path']; assert digest(ip) == cfg['source_identity_review']['sha256']
    identity = json.loads(ip.read_text(encoding='utf-8')); crosswalk = identity['source_origins']
    for origin, source in crosswalk.items():
        chosen = obj.obs.origin.astype(str).eq(origin)
        assert chosen.any() and obj.obs.loc[chosen, 'donor'].astype(str).eq(source['source_author_donor']).all()
        assert obj.obs.loc[chosen, 'Meta:Vaccination_Status'].astype(str).eq(source['source_phase']).all()
    obj.obs['donor_id'] = obj.obs.donor.copy()
    obj.obs['source_GEO_accession'] = pd.Categorical([crosswalk[o]['source_GSM'] for o in obj.obs.origin.astype(str)])
    obj.obs['source_GEO_individual_as_reported'] = pd.Categorical([crosswalk[o]['original_GEO_individual'] for o in obj.obs.origin.astype(str)])
    obj.obs['source_condition_as_reported'] = pd.Categorical([crosswalk[o]['source_cohort_group'] for o in obj.obs.origin.astype(str)])
    obj.obs['source_visit_as_reported'] = obj.obs['Meta:Vaccination_Status'].copy()
    obj.obs['tissue'] = 'blood'; obj.obs['is_primary_data'] = True
    obj.obs['source_scope'] = audit['source_scope']
    obj.obs['source_identity_policy'] = audit['identity_policy']
    obj.obs['source_time_policy'] = 'Literal pre/post vaccination and post-infection source labels; no inferred day, dose or calendar clock.'
    obj.obs['source_doublet_policy'] = 'Annotation-table membership is retained as provenance; no biological doublet classification inferred from identity columns.'
    obj.obs['source_split_assignment'] = 'Unassigned; all visits of one source participant stay together. Global cross-study aliases unresolved.'
    for column in obj.obs:
        if obj.obs[column].dtype == object: obj.obs[column] = pd.Categorical(obj.obs[column])
    needed = obj.X.data.nbytes + obj.X.indices.nbytes + obj.X.indptr.nbytes
    if footprint() + 2 * needed + 2_000_000_000 >= LIMIT: raise RuntimeError('Native RNA staging exceeds project disk budget')
    obj.write_h5ad(destination, compression='gzip')
    return dict(kind='original_HSCT_PBMC_exact_author_identity_subset', native_count_audit=cfg['native_count_audit'],
        source_identity_review=cfg['source_identity_review'], source_cells=136901, explicitly_identified_cells=89753,
        unassigned_source_RNA_cells=47148, source_origins=16, source_participant_codes=7,
        arithmetic='Unchanged original RNA rows selected by literal author barcode identity; no inferred suffix mapping, normalization, gene intersection or averaging.'), [ROOT / a['path'] for a in audit['original_assets'].values()] + [bundle]
