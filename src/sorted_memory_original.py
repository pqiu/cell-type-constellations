"""Stage unchanged sorted PBMC RNA with exact source cell and gene identities."""
import json

import anndata as ad
import pandas as pd

from scripts.download_pbmc_asset import ROOT, LIMIT, footprint
from src.rag_pbmc_original import verify_arrays


def prepare(cfg, destination):
    from src.pbmc_pipeline import digest
    from scripts.prepare_PBMC_106540 import replay_sources
    audit = replay_sources()
    assert digest(ROOT / cfg['native_count_audit']['path']) == cfg['native_count_audit']['sha256']
    assert audit['complete'] and audit['all_native_versioned_gene_IDs_exactly_match_GENCODE19']
    bundle = ROOT / audit['native_bundle']['path']
    if not bundle.exists():
        # Original validation removed the staged count carrier. Recreate it
        # from the pinned raw files, then put provenance records back byte-for-byte.
        from scripts.audit_PBMC_106540_original_counts import main as rebuild
        code = ROOT / 'scripts/audit_PBMC_106540_original_counts.py'
        assert digest(code) == audit['code_sha256']
        preserved = [ROOT / cfg['native_count_audit']['path'],
                     ROOT / audit['original_GEO_cell_records']['path'],
                     ROOT / 'references/datasets/GSE106540_original_source_review_manifest.json']
        old_bytes = {path: path.read_bytes() for path in preserved if path.exists()}
        try:
            regenerated = rebuild()
            assert {k: v for k, v in regenerated.items() if k != 'checked_at'} == {
                k: v for k, v in audit.items() if k != 'checked_at'}
        finally:
            for path, raw in old_bytes.items():
                path.write_bytes(raw)
    assert digest(bundle) == audit['native_bundle']['sha256']
    obj = ad.read_h5ad(bundle); verify_arrays(obj.X, audit['native_CSR_array_pins'])
    assert obj.shape == (2244, 57820) and obj.obs_names.is_unique
    ip = ROOT / cfg['source_identity_review']['path']; assert digest(ip) == cfg['source_identity_review']['sha256']
    identity = json.loads(ip.read_text(encoding='utf-8'))
    for view, source in identity['source_views'].items():
        chosen = obj.obs.donor_id.astype(str).eq(source['donor_id']) & obj.obs.source_sorted_subtype.astype(str).eq(source['source_sorted_subtype'])
        assert chosen.sum() == source['original_cells']
        assert sorted(obj.obs.loc[chosen, 'source_GSM'].astype(str)) == source['original_GSMs']
    obj.obs['tissue'] = 'blood'; obj.obs['is_primary_data'] = True
    obj.obs['source_scope'] = audit['source_scope']
    obj.obs['source_identity_policy'] = audit['identity_policy']
    obj.obs['source_history_policy'] = audit['history_policy']
    obj.obs['source_time_policy'] = 'Exact visit/aliquot relationship not reported; sort populations remain dependent views of one conservative source participant/material group.'
    obj.obs['source_split_assignment'] = 'Unassigned; all source views from one participant stay together. Global cross-study aliases unresolved.'
    for col in obj.obs:
        if obj.obs[col].dtype == object: obj.obs[col] = pd.Categorical(obj.obs[col])
    needed = obj.X.data.nbytes + obj.X.indices.nbytes + obj.X.indptr.nbytes
    if footprint() + 2 * needed + 2_000_000_000 >= LIMIT: raise RuntimeError('Native RNA staging exceeds project disk budget')
    obj.write_h5ad(destination, compression='gzip')
    return dict(kind='original_sorted_memory_PBMC_RNA', native_count_audit=cfg['native_count_audit'],
        source_identity_review=cfg['source_identity_review'], source_cells=2244, native_genes=57820,
        source_participant_material_groups=16, dependent_sorted_views=28,
        arithmetic='Unchanged original HTSeq-read counts, versioned gene IDs and exact GENCODE19symbols; no normalized values or averaging used during staging.'), [ROOT / a['path'] for a in audit['original_assets'].values()] + [bundle]
