"""Stage native memory-B RNA without imposing the smaller publication selection."""
import json

import anndata as ad
import pandas as pd

from scripts.download_pbmc_asset import ROOT, LIMIT, footprint
from src.rag_pbmc_original import verify_arrays


def prepare(cfg,destination):
    from src.pbmc_pipeline import digest
    from scripts.prepare_PBMC_157447 import replay_sources
    audit=replay_sources();assert digest(ROOT/cfg['native_count_audit']['path'])==cfg['native_count_audit']['sha256']
    assert audit['complete'] and audit['all_original_stored_values_finite_nonnegative_integer']
    bundle=ROOT/audit['native_bundle']['path']
    if not bundle.exists():
        from scripts.audit_PBMC_157447_original_counts import main as rebuild
        code=ROOT/'scripts/audit_PBMC_157447_original_counts.py'
        assert digest(code)==audit['code_sha256']
        audit_path=ROOT/cfg['native_count_audit']['path']
        manifest_paths=[ROOT/'references/datasets/GSE157447_resolved_original_source_manifest.json',
                        ROOT/'references/datasets/GSE157447_original_source_review_manifest.json']
        old_bytes={path:path.read_bytes() for path in [audit_path,*manifest_paths] if path.exists()}
        try:
            regenerated=rebuild()
            assert {k:v for k,v in regenerated.items() if k!='checked_at'}=={
                k:v for k,v in audit.items() if k!='checked_at'}
        finally:
            for path,raw in old_bytes.items():path.write_bytes(raw)
    assert digest(bundle)==audit['native_bundle']['sha256']
    obj=ad.read_h5ad(bundle);verify_arrays(obj.X,audit['native_CSR_array_pins']);assert obj.shape==(825,58051) and obj.obs_names.is_unique
    ip=ROOT/cfg['source_identity_review']['path'];assert digest(ip)==cfg['source_identity_review']['sha256']
    identity=json.loads(ip.read_text(encoding='utf-8'))
    for view,source in identity['source_views'].items():
        chosen=obj.obs.donor_id.astype(str).eq(source['donor_id']) & obj.obs.source_plate_id.astype(str).eq(source['source_plate_id']) & obj.obs.source_title_time_token.astype(str).eq(source['source_title_time_token'])
        assert chosen.sum()==source['original_cells'] and sorted(obj.obs.loc[chosen,'source_GSM'].astype(str))==source['original_GSMs']
    obj.obs['tissue']='blood';obj.obs['is_primary_data']=True
    obj.obs['source_scope']=audit['source_scope'];obj.obs['source_identity_policy']=audit['identity_policy']
    obj.obs['source_time_policy']=audit['time_policy'];obj.obs['source_author_selection_policy']=identity['author_selection_policy']
    obj.obs['source_split_assignment']='Unassigned; all source views from one donor stay together. Global cross-study aliases unresolved.'
    for col in obj.obs:
        if obj.obs[col].dtype==object:obj.obs[col]=pd.Categorical(obj.obs[col])
    needed=obj.X.data.nbytes+obj.X.indices.nbytes+obj.X.indptr.nbytes
    if footprint()+2*needed+2_000_000_000>=LIMIT:raise RuntimeError('Native RNA staging exceeds project disk budget')
    obj.write_h5ad(destination,compression='gzip')
    return dict(kind='original_memory_B_PBMC_read_counts',native_count_audit=cfg['native_count_audit'],source_identity_review=cfg['source_identity_review'],
        original_cells=825,native_genes=58051,source_participant_material_groups=12,dependent_source_views=20,
        arithmetic='Unchanged original STAR read counts and native symbols; all provided source cells before ownQC, without imposing source publication plate/gene selection.'),[ROOT/a['path'] for a in audit['original_assets'].values()]+[bundle]
