"""Stage one complete native gene group from reviewed original PBMC counts."""
import json

import anndata as ad

from scripts.download_pbmc_asset import ROOT, LIMIT, footprint
from src.rag_pbmc_original import verify_arrays


def prepare(cfg, destination):
    from scripts.prepare_PBMC_180578 import replay_sources, checked, checked_native_bundle
    audit = replay_sources(cfg['dataset_id'])
    checked(cfg['native_count_audit'])
    identity = json.loads(checked(cfg['source_identity_review']).read_text(encoding='utf-8'))
    assert audit['source_specimen_identity_review'] == cfg['source_identity_review']
    assert identity['source_identity_review_complete']
    assert audit['complete'] and audit['all_original_stored_values_finite_nonnegative_integer']
    bundle = checked_native_bundle(audit['native_bundle'])
    obj = ad.read_h5ad(bundle)
    verify_arrays(obj.X, audit['native_CSR_array_pins'])
    assert obj.shape == (audit['original_cells'], audit['original_genes'])
    assert obj.obs_names.is_unique and obj.var_names.is_unique
    assert obj.obs.tissue.astype(str).eq('blood').all() and obj.obs.is_primary_data.eq(True).all()
    assert obj.obs.source_native_gene_coordinate_group.astype(str).eq(cfg['dataset_id']).all()
    assert obj.obs.source_specimen_id.nunique() == audit['source_draws']
    needed = obj.X.data.nbytes + obj.X.indices.nbytes + obj.X.indptr.nbytes
    if footprint() + 2 * needed + 2_000_000_000 >= LIMIT:
        raise RuntimeError('Original RNA staging exceeds project disk budget')
    obj.write_h5ad(destination, compression='gzip')
    return dict(kind='reviewed_original_PBMC_native_gene_group', native_count_audit=cfg['native_count_audit'],
                source_identity_review=cfg['source_identity_review'], original_cells=obj.n_obs,
                native_genes=obj.n_vars, arithmetic='Unchanged original counts, ordered axes and source metadata.'), [bundle]
