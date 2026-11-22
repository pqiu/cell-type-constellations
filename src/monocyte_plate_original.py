"""Stage unchanged native RNA with exact plate-to-GEO source identities."""
import json

import anndata as ad
import pandas as pd

from scripts.download_pbmc_asset import ROOT, LIMIT, footprint
from src.rag_pbmc_original import verify_arrays


def prepare(cfg, destination):
    from src.pbmc_pipeline import digest
    from scripts.prepare_PBMC_184241 import replay_sources
    audit = replay_sources()
    assert digest(ROOT / cfg['native_count_audit']['path']) == cfg['native_count_audit']['sha256']
    assert audit['complete'] and audit['all_original_stored_values_finite_nonnegative_integer']
    bundle = ROOT / audit['native_bundle']['path']; assert bundle.exists()
    obj = ad.read_h5ad(bundle); verify_arrays(obj.X, audit['native_CSR_array_pins'])
    assert obj.shape == (1710, 23895) and obj.obs_names.is_unique
    ip = ROOT / cfg['source_identity_review']['path']; assert digest(ip) == cfg['source_identity_review']['sha256']
    identity = json.loads(ip.read_text(encoding='utf-8'))
    for gsm, source in identity['source_views'].items():
        chosen = obj.obs.source_GSM.astype(str).eq(gsm); assert chosen.sum() == 95
        assert obj.obs.loc[chosen, 'donor_id'].astype(str).eq(source['source_donor']).all()
        assert obj.obs.loc[chosen, 'source_visit_as_reported'].astype(str).eq(source['source_time']).all()
        assert obj.obs.loc[chosen, 'source_treatment_as_reported'].astype(str).eq(source['source_treatment']).all()
    obj.obs['tissue'] = 'blood'; obj.obs['is_primary_data'] = True
    obj.obs['source_condition_as_reported'] = 'Healthy donor (source series design)'
    obj.obs['source_scope'] = audit['source_scope']
    obj.obs['source_identity_policy'] = audit['sample_identity_policy']
    obj.obs['source_time_policy'] = 'Preserve literal original visit labels; ex vivo condition is distinct from clinical condition.'
    obj.obs['source_split_assignment'] = 'Unassigned; all visits and ex vivo views from one donor stay together; cross-study aliases unresolved.'
    for column in obj.obs:
        if obj.obs[column].dtype == object: obj.obs[column] = pd.Categorical(obj.obs[column])
    needed = obj.X.data.nbytes + obj.X.indices.nbytes + obj.X.indptr.nbytes
    if footprint() + 2 * needed + 2_000_000_000 >= LIMIT: raise RuntimeError('Native RNA staging exceeds project disk budget')
    obj.write_h5ad(destination, compression='gzip')
    return dict(kind='original_monocyte_plate_RNA', native_count_audit=cfg['native_count_audit'],
        source_identity_review=cfg['source_identity_review'], source_cells=1710, native_genes=23895,
        source_visits=9, participant_codes=3, dependent_views=18, mitochondrial_QC_available=False,
        arithmetic='Unchanged native RNA rows and symbols; exact cell and plate identities, no normalization or averaging during staging.'), [ROOT / a['path'] for a in audit['original_assets'].values()] + [bundle]
