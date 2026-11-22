"""Stage a pinned public native-count bundle without assuming a GEO accession."""
import importlib
import hashlib
import json
import re
import shutil
import anndata as ad
from scripts.download_pbmc_asset import ROOT, LIMIT, footprint
from src.rag_pbmc_original import verify_arrays


def checked(pin):
    from src.pbmc_pipeline import digest
    path = (ROOT / pin['path']).resolve()
    assert path.is_relative_to(ROOT.resolve()) and path.is_file()
    assert digest(path) == pin['sha256']
    if 'bytes' in pin: assert path.stat().st_size == pin['bytes']
    return path


def replay_module(name):
    if not re.fullmatch(r'scripts\.prepare_PBMC_[A-Za-z0-9_]+', name):
        raise ValueError('Expected a local reviewed PBMC replay module')
    return name.replace('.', '/') + '.py'


def stage_bundle(cfg, audit, destination):
    assert audit['complete'] and audit['all_original_stored_values_finite_nonnegative_integer']
    identity = json.loads(checked(cfg['source_identity_review']).read_text(encoding='utf-8'))
    assert identity['original_count_audit'] == cfg['native_count_audit']
    bundle = checked(audit['native_bundle'])
    destination = destination.resolve()
    assert destination.is_relative_to((ROOT/'data/downloads').resolve()) and destination != bundle
    obj = ad.read_h5ad(bundle)
    verify_arrays(obj.X, audit['native_CSR_array_pins'])
    assert obj.shape == (audit['original_cells'],audit['original_genes'])
    assert obj.obs_names.is_unique and obj.var_names.is_unique
    assert obj.obs.tissue.astype(str).eq('blood').all() and obj.obs.is_primary_data.eq(True).all()
    if destination.exists():
        # A previous pinned replay may have left this identical native carrier.
        # Reuse only exact bytes; never overwrite or silently accept a variant.
        with bundle.open('rb') as source, destination.open('rb') as staged:
            if hashlib.file_digest(source,'sha256').digest() != hashlib.file_digest(staged,'sha256').digest():
                raise ValueError('Existing native staging differs from the pinned source bundle')
    else:
        needed = obj.X.data.nbytes + obj.X.indices.nbytes + obj.X.indptr.nbytes
        if footprint() + 2*needed + 2_000_000_000 >= LIMIT or shutil.disk_usage(ROOT).free < 2*needed + 2_000_000_000:
            raise RuntimeError('Native RNA staging exceeds available project disk budget')
        temporary = destination.with_suffix('.building.h5ad')
        assert not temporary.exists(), 'Inspect partial staging before retrying'
        obj.write_h5ad(temporary,compression='gzip');temporary.replace(destination)
    manifest = dict(kind='reviewed_public_native_RNA_bundle', native_count_audit=cfg['native_count_audit'],
                    source_identity_review=cfg['source_identity_review'], original_cells=obj.n_obs,native_genes=obj.n_vars,
                    arithmetic='Unchanged original RNA counts and explicit native gene/cell/source axes. No normalization, intersection, missing-gene zero filling or averaging during staging.')
    return manifest, [checked(p) for p in cfg['reviewed_source_assets']] + [bundle]


def prepare(cfg,destination):
    path = replay_module(cfg['native_replay_module'])
    assert cfg['native_replay_code']['path'] == path
    checked(cfg['native_replay_code'])
    module=importlib.import_module(cfg['native_replay_module'])
    if cfg['dataset_id']=='GSE189050':
        audit_path=ROOT/cfg['native_count_audit']['path']
        pinned=json.loads(audit_path.read_text(encoding='utf-8'))
        bundle=ROOT/pinned['native_bundle']['path']
        if not bundle.exists():
            from scripts.rebuild_gse189050_native_bundle import main as rebuild_guarded
            rebuild_guarded()
    if cfg['dataset_id'].startswith('GSE329190_'):
        from scripts.replay_PBMC_329190_partition_cache import main as restore_partition_cache
        restore_partition_cache(cfg['dataset_id'].split('_',1)[1])
    if cfg['dataset_id'] in ('GSE290722','GSE319710','PRJEB86440') or cfg['dataset_id'].startswith('GSE329190_'):
        pinned=json.loads(checked(cfg['native_count_audit']).read_text(encoding='utf-8'))
        bundle=ROOT/pinned['native_bundle']['path']
        run=ROOT/'results/star_maps'/cfg['dataset_id']/'run.json'
        held=run.with_name('run.hierarchy-replay-hold.json')
        if not bundle.exists() and run.exists():
            assert not held.exists()
            run.rename(held)
            try:audit=module.replay_sources()
            finally:
                if held.exists():held.rename(run)
        else:audit=module.replay_sources()
    else:audit=module.replay_sources()
    assert json.loads(checked(cfg['native_count_audit']).read_text(encoding='utf-8')) == audit
    return stage_bundle(cfg,audit,destination)
