"""Stage an immutable, reviewed original RNA bundle without expression changes."""
import importlib
import json
import re

import anndata as ad

from scripts.download_pbmc_asset import ROOT,LIMIT,footprint
from src.rag_pbmc_original import verify_arrays


def prepare(cfg,destination):
    from src.pbmc_pipeline import digest
    study=cfg['study_id'];assert re.fullmatch(r'GSE[0-9]+',study)
    replay=importlib.import_module('scripts.prepare_PBMC_'+study[3:]).replay_sources
    if study in ('GSE157278','GSE167976'):
        # Its historical replay insists on byte identity of a cleaned H5AD
        # container. The reconstructed carrier has independently pinned arrays
        # and axes but different container bytes; audit source assets directly.
        from scripts.download_pbmc_asset import download
        audit=json.loads((ROOT/cfg['native_count_audit']['path']).read_text(encoding='utf-8'))
        for asset in audit['original_assets'].values():
            path=download(asset['url'],ROOT/asset['path'],asset['bytes'])
            assert digest(path)==asset['sha256']
    elif study in ('GSE305477','GSE336230','GSE345981'):
        pinned=json.loads((ROOT/cfg['native_count_audit']['path']).read_text(encoding='utf-8'))
        bundle=ROOT/pinned['native_bundle']['path']
        run=ROOT/'results/star_maps'/study/'run.json'
        held=run.with_name('run.hierarchy-replay-hold.json')
        if not bundle.exists() and run.exists():
            assert not held.exists()
            run.rename(held)
            try:audit=replay()
            finally:
                if held.exists():held.rename(run)
        else:audit=replay()
    elif study in ('GSE295173','GSE295174'):
        # Historical overlap audits scan completed run records by literal text.
        # Two unrelated curated configs contain copied GSE295174 prose and
        # this source's own run also self-matches. Preserve the true companion
        # GSE295173 match when rebuilding GSE295174.
        audit_path=ROOT/cfg['native_count_audit']['path']
        pinned=json.loads(audit_path.read_text(encoding='utf-8'))
        bundle=ROOT/pinned['native_bundle']['path']
        if bundle.exists():
            audit=replay()
        else:
            hidden=['GSE276742','GSE324824',study]
            if study=='GSE295173':hidden.append('GSE295174')
            moved=[]
            try:
                for other in hidden:
                    run=ROOT/'results/star_maps'/other/'run.json'
                    held=run.with_name('run.hierarchy-replay-hold.json')
                    assert not held.exists()
                    if run.exists():run.rename(held);moved.append((held,run))
                audit=replay()
            finally:
                for held,run in reversed(moved):held.rename(run)
    else:
        audit=replay()
    assert digest(ROOT/cfg['native_count_audit']['path'])==cfg['native_count_audit']['sha256']
    assert audit['complete'] and audit['all_original_stored_values_finite_nonnegative_integer']
    ip=ROOT/cfg['source_identity_review']['path'];assert digest(ip)==cfg['source_identity_review']['sha256']
    identity=json.loads(ip.read_text(encoding='utf-8'));assert identity['original_count_audit']==cfg['native_count_audit']
    bundle=ROOT/audit['native_bundle']['path']
    if not bundle.exists():
        if study=='GSE157278':
            from scripts.rebuild_gse157278_native_bundle import main as rebuild_guarded
            rebuild_guarded()
        elif study=='GSE167976':
            from scripts.rebuild_gse167976_native_bundle import main as rebuild_guarded
            rebuild_guarded()
        elif study=='GSE175450':
            # This source pins its staging code as an asset rather than a
            # top-level code_sha256 and uses a reviewed sparse-count cache.
            from scripts.prepare_PBMC_180578 import checked
            from scripts.audit_PBMC_175450_original_counts import main as rebuild_counts
            from scripts.stage_PBMC_175450_native_counts import main as rebuild_bundle
            checked(audit['code'])
            review=audit['original_count_review']
            if not (ROOT/audit['derived_author_assigned_native_RNA_cache']['path']).exists():
                assert digest(ROOT/'scripts/audit_PBMC_175450_original_counts.py')==json.loads((ROOT/review['path']).read_text())['code']['sha256']
                rebuilt_review=rebuild_counts()
                assert digest(ROOT/review['path'])==review['sha256']
            regenerated=rebuild_bundle()
            assert regenerated==audit
            assert digest(ROOT/cfg['native_count_audit']['path'])==cfg['native_count_audit']['sha256']
        elif study=='GSE199994':
            from scripts.rebuild_gse199994_native_bundle import main as rebuild_guarded
            rebuild_guarded()
        else:
            # The original validation cleaned this staged carrier. Rebuild from
            # source assets and retain the original provenance record byte-for-byte.
            if study=='GSE158769':
                from scripts.replay_PBMC_158769_sparse_blocks import main as replay_blocks
                replay_blocks()
            if study=='GSE235857':
                from scripts.replay_PBMC_235857_sparse_materials import main as replay_materials
                replay_materials()
            audit_code=(ROOT/'scripts/prepare_PBMC_248256.py' if study=='GSE248256'
                        else ROOT/f'scripts/audit_PBMC_{study[3:]}_original_counts.py')
            assert digest(audit_code)==audit['code_sha256']
            rebuild=importlib.import_module(f'scripts.audit_PBMC_{study[3:]}_original_counts').main
            audit_path=ROOT/cfg['native_count_audit']['path']
            manifest_paths=[ROOT/f'references/datasets/{study}_resolved_original_source_manifest.json',
                            ROOT/f'references/datasets/{study}_original_source_review_manifest.json']
            old_bytes={path:path.read_bytes() for path in [audit_path,*manifest_paths] if path.exists()}
            try:
                regenerated=rebuild()
                assert {k:v for k,v in regenerated.items() if k!='checked_at'}=={
                    k:v for k,v in audit.items() if k!='checked_at'}
            finally:
                for path,raw in old_bytes.items():path.write_bytes(raw)
    if digest(bundle)!=audit['native_bundle']['sha256']:
        if study not in ('GSE157278','GSE167976','GSE199994'):raise ValueError('Rebuilt native carrier differs from original audited bytes')
        proof_path=ROOT/'results/hierarchy_star_maps'/study/'native_rebuild_proof.json'
        proof=json.loads(proof_path.read_text(encoding='utf-8'))
        assert proof['audit_sha256']==digest(ROOT/cfg['native_count_audit']['path'])
        assert proof['source_cell_metadata_sha256']==audit['source_cell_metadata']['sha256']
        assert proof['native_gene_coordinates_sha256']==audit['native_gene_coordinates']['sha256']
        assert proof['native_CSR_array_pins']==audit['native_CSR_array_pins']
        assert proof['rebuilt_h5ad_sha256']==digest(bundle)
    obj=ad.read_h5ad(bundle);verify_arrays(obj.X,audit['native_CSR_array_pins'])
    assert obj.shape==(audit['original_cells'],audit['original_genes']) and obj.obs_names.is_unique and obj.var_names.is_unique
    assert obj.obs.tissue.astype(str).eq('blood').all() and obj.obs.is_primary_data.eq(True).all()
    needed=obj.X.data.nbytes+obj.X.indices.nbytes+obj.X.indptr.nbytes
    if footprint()+2*needed+2_000_000_000>=LIMIT:raise RuntimeError('Native RNA staging exceeds project disk budget')
    obj.write_h5ad(destination,compression='gzip')
    return dict(kind='reviewed_original_native_RNA_bundle',native_count_audit=cfg['native_count_audit'],source_identity_review=cfg['source_identity_review'],
        original_cells=audit['original_cells'],native_genes=audit['original_genes'],
        arithmetic='Unchanged reviewed original native RNA counts, axes and source metadata. No normalization, gene intersection or averaging during staging.'),[ROOT/a['path'] for a in audit['original_assets'].values()]+[bundle]
