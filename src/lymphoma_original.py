"""Preserve complete original RNA and author metadata before reviewed selection."""
import json

import anndata as ad

from scripts.download_pbmc_asset import ROOT,LIMIT,footprint
from scripts.prepare_PBMC_180578 import checked
from src.rag_pbmc_original import verify_arrays


def prepare(cfg,destination):
    from scripts.prepare_PBMC_207192 import replay_sources
    audit=replay_sources();checked(cfg['native_count_audit'])
    identity=json.loads(checked(cfg['source_identity_review']).read_text(encoding='utf-8'))
    assert identity['original_count_audit']==cfg['native_count_audit'] and identity['eligible_for_sample_checkpoint']
    assert audit['complete_original_post_RNA_value_and_author_identity_audit'] and audit['all_original_stored_values_finite_nonnegative_integer']
    obj=ad.read_h5ad(checked(audit['native_bundle']));verify_arrays(obj.X,audit['native_CSR_array_pins'])
    assert obj.shape==(audit['original_cells'],audit['original_genes']) and obj.obs_names.is_unique and obj.var_names.is_unique
    obj.obs['source_identity_status']=['held unresolved capture provenance' if c in identity['excluded_unresolved_captures'] else 'reviewed original author participant/time' for c in obj.obs.source_capture_label.astype(str)]
    obj.obs['source_split_assignment']='Unassigned; keep all visits for one source participant together; cross-study aliases unresolved.'
    needed=obj.X.data.nbytes+obj.X.indices.nbytes+obj.X.indptr.nbytes
    if footprint()+2*needed+2_000_000_000>=LIMIT:raise RuntimeError('Original RNA stage exceeds project disk budget')
    obj.write_h5ad(destination,compression='gzip')
    return dict(kind='original_post_infusion_native_RNA',native_count_audit=cfg['native_count_audit'],source_identity_review=cfg['source_identity_review'],original_cells=obj.n_obs,native_genes=obj.n_vars,arithmetic='Unchanged native RNA counts and axes; original author metadata retained; reviewed selection occurs separately.'),[ROOT/a['path'] for a in audit['original_assets'].values()]+[ROOT/audit['native_bundle']['path']]
