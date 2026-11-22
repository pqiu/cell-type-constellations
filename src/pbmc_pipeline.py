"""Conservative, blockwise CELLxGENE raw-count annotation and star export.

Dataset configurations must explicitly identify biological specimens. No automatic
donor-to-sample substitution or inferred PBMC restriction is permitted.
"""
import hashlib
import importlib.metadata
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.download_pbmc_asset import download, footprint, LIMIT

REPRESENTATIONS = ('raw', 'norm10k', 'norm10k_log1p')
DEFAULT_OUTPUT_KEYS = tuple(
    rep + '_' + stat
    for rep in REPRESENTATIONS
    for stat in ('mean', 'variance', 'sem')
) + ('raw_sum',)


def output_statistic_keys(cfg):
    """Validate an explicit subset of persisted star statistics."""
    requested = cfg.get('star_output_statistics')
    if requested is None:
        return DEFAULT_OUTPUT_KEYS
    if not isinstance(requested, list) or not requested or len(requested) != len(set(requested)):
        raise ValueError('star_output_statistics must be a nonempty unique list')
    unknown = [key for key in requested if key not in DEFAULT_OUTPUT_KEYS]
    if unknown:
        raise ValueError('Unknown star output statistics: ' + ','.join(unknown))
    if 'norm10k_log1p_mean' not in requested:
        raise ValueError('Persisted star geometry requires norm10k_log1p_mean')
    reason = cfg.get('star_output_statistics_reason')
    if not isinstance(reason, str) or not reason.strip():
        raise ValueError('A reduced statistic set requires an explicit reason')
    return tuple(requested)

def donor_identity_column(cfg):
    """Return the reviewed specimen-donor field without silently reinterpreting donor_id."""
    column = cfg.get('donor_identity_column', 'donor_id')
    if not isinstance(column, str) or not column.strip():
        raise ValueError('donor_identity_column must be a nonempty source column name')
    if column != 'donor_id':
        evidence = cfg.get('donor_identity_mapping_evidence')
        if not isinstance(evidence, str) or not evidence.strip():
            raise ValueError('A non-default donor identity column requires reviewed mapping evidence')
    return column

def summarize(matrix):
    """Unbiased cell variance and conditional SEM; n=1 uncertainty is undefined."""
    x = sparse.csr_matrix(matrix, dtype=np.float64)
    n = x.shape[0]
    if n == 0:
        raise ValueError('Cannot summarize an absent type')
    total = np.asarray(x.sum(axis=0)).ravel()
    mean = total / n
    if n > 1:
        square = np.asarray(x.multiply(x).sum(axis=0)).ravel()
        var = np.maximum((square - total * mean) / (n - 1), 0)
        sem = np.sqrt(var / n)
    else:
        var = sem = np.full(x.shape[1], np.nan)
    return mean, var, sem, total

def representations(counts):
    x = sparse.csr_matrix(counts, dtype=np.float64)
    totals = np.asarray(x.sum(axis=1)).ravel()
    if np.any(totals <= 0):
        raise ValueError('Zero-library cells must be excluded before normalization')
    norm = x.multiply((10000 / totals)[:, None]).tocsr()
    log = norm.copy(); log.data = np.log1p(log.data)
    return x, norm, log

def gene_alignment(symbols, gene_order):
    target = {g: i for i, g in enumerate(gene_order)}
    pairs = [(i, target[g]) for i, g in enumerate(symbols) if g in target]
    if not pairs:
        raise ValueError('No genes overlap the annotation reference')
    rows, cols = zip(*pairs)
    return sparse.csr_matrix((np.ones(len(pairs)), (rows, cols)), shape=(len(symbols), len(gene_order)))

def annotation_symbols(var, cfg):
    """Return configured gene symbols, with an optional explicit source prefix removed."""
    symbols = (var.index.astype(str).to_numpy() if cfg['annotation_gene_column'] == '__index__'
               else var[cfg['annotation_gene_column']].astype(str).to_numpy())
    prefix = cfg.get('annotation_gene_prefix_to_strip')
    if prefix is None:
        return symbols
    if not isinstance(prefix, str) or not prefix:
        raise ValueError('annotation_gene_prefix_to_strip must be a nonempty string')
    source = symbols.astype(str)
    if not any(symbol.startswith(prefix) for symbol in source):
        raise ValueError('Configured annotation gene prefix is absent from the source gene axis')
    # Mixed host/pathogen references are valid: strip the host-build prefix where
    # present and retain explicitly named non-host features for provenance/QC.
    return np.asarray([symbol.removeprefix(prefix) for symbol in source])

def annotation_input(counts, mapper):
    """Match SCimilarity: align raw counts, normalize aligned genes, then log1p."""
    aligned = (sparse.csr_matrix(counts, dtype=np.float64) @ mapper).tocsr()
    library = np.asarray(aligned.sum(axis=1)).ravel()
    if np.any(library <= 0):
        raise ValueError('A cell has no counts in the annotation gene space')
    normalized = aligned.multiply((10000 / library)[:, None]).tocsr()
    normalized.data = np.log1p(normalized.data)
    return normalized, library

def digest(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()

def validate_outputs(out):
    report = json.loads((out / 'run.json').read_text(encoding='utf-8'))
    for name, checksum in report['output_sha256'].items():
        if digest(out / name) != checksum:
            raise ValueError('Output checksum mismatch: ' + name)
    stars = pd.read_csv(out / 'stars.csv')
    genes = pd.read_csv(out / 'genes.csv')
    # Source identity components are nominal strings; CSV inference must not
    # turn a source participant "01" into "1" during identity validation.
    from src.specimen_views import read_samples_preserving_identity
    samples = read_samples_preserving_identity(out / 'samples.csv', report['configuration'])
    if report['configuration'].get('physical_sample_columns') is not None:
        from src.specimen_views import attach_physical_identity
        if 'physical_specimen_uid' not in samples:raise ValueError('Missing physical specimen identity')
        expected=attach_physical_identity(samples.drop(columns=['physical_specimen_uid']),report['configuration'])
        if not samples.physical_specimen_uid.eq(expected.physical_specimen_uid).all():
            raise ValueError('Physical specimen identity differs from source columns')
        if report['biological_samples']!=samples.physical_specimen_uid.nunique() or report.get('observation_views')!=samples.sample_id.nunique():
            raise ValueError('Physical specimen/view counts do not reconcile')
    title = report['configuration'].get('source_metadata', {}).get('dataset_metadata', {}).get('title', '')
    day = re.search(r'(\d+) days? post-intubation', title)
    if day and 'days_post_intubation' in samples and not samples.days_post_intubation.astype(str).eq(day.group(1)).all():
        raise ValueError('Source title and cell-level collection time disagree; quarantine until resolved')
    z = h5py.File(out / 'expression.h5', 'r') if (out / 'expression.h5').exists() else np.load(out / 'expression.npz', allow_pickle=False)
    output_keys = output_statistic_keys(report['configuration'])
    optional_support = {'gene_n_cells','source_gene_availability','source_library_names'}
    if set(z.keys()) - optional_support != set(output_keys):
        raise ValueError('Persisted expression statistics differ from configuration')
    if report['configuration'].get('use_source_gene_availability', False):
        required = ['gene_n_cells','source_gene_availability','source_library_names']
        if any(name not in z for name in required):
            raise ValueError('Missing source gene availability output')
        availability = np.asarray(z['source_gene_availability'])
        names = np.asarray(z['source_library_names'])
        if availability.dtype != np.dtype(bool) or availability.shape != (len(names),len(genes)) or len(set(names)) != len(names):
            raise ValueError('Invalid saved source-library gene availability')
        if z['gene_n_cells'].shape != (len(stars),len(genes)):
            raise ValueError('Invalid gene support dimensions')
        support_dataset = z['gene_n_cells']
        if not isinstance(support_dataset,h5py.Dataset) or support_dataset.chunks != (1,len(genes)) or support_dataset.fillvalue != 0:
            raise ValueError('Invalid gene support storage convention')
        allocated_support = {support_dataset.id.get_chunk_info(i).chunk_offset[0] for i in range(support_dataset.id.get_num_chunks())}
        integer_availability = availability.astype(np.int64)
        cells = pd.read_csv(out/'cell_metadata.csv.gz',usecols=['retained','sample_id','scimilarity_label','source_library_index'])
        cells = cells.loc[cells.retained.eq(True)]
        libraries = cells.source_library_index.to_numpy()
        if libraries.dtype.kind not in 'iu' or (libraries < 0).any() or (libraries >= len(names)).any():
            raise ValueError('Invalid retained source-library identities')
        grouped = cells.groupby(['sample_id','scimilarity_label','source_library_index']).size()
        frequencies = {}
        for (sid,label,library),n in grouped.items():
            frequencies.setdefault((sid,label),np.zeros(len(names),dtype=np.int64))[library] = n
        for i,star in stars.iterrows():
            stored_frequencies = frequencies.get((star.sample_id,star.label))
            if star.n_cells == 0 and i not in allocated_support and stored_frequencies is None:
                continue
            weights = stored_frequencies if stored_frequencies is not None else np.zeros(len(names),dtype=np.int64)
            if weights.sum() != star.n_cells or not np.array_equal(z['gene_n_cells'][i],weights @ integer_availability):
                raise ValueError('Invalid per-gene observed-cell support')
    for rep in REPRESENTATIONS:
        if rep + '_mean' not in output_keys:
            continue
        a = z[rep + '_mean']
        if a.shape != (len(stars), len(genes)):
            raise ValueError('Star/gene dimensions do not match metadata')
        supported = stars.status.eq('supported').to_numpy()
        if isinstance(a, h5py.Dataset):
            if a.chunks != (1, len(genes)) or not np.isnan(a.fillvalue):
                raise ValueError('Unexpected chunk or missingness convention')
            allocated = {a.id.get_chunk_info(i).chunk_offset[0] for i in range(a.id.get_num_chunks())}
            # Unallocated rows are guaranteed NaN by the HDF5 fill value. Inspect
            # every allocated row and require storage for every supported star.
            if not set(np.flatnonzero(supported)).issubset(allocated):
                raise ValueError('A supported star has no stored centroid')
            indices = sorted(allocated)
        else:
            indices = range(len(stars))
        for i in indices:
            values = a[i]
            expected_finite = np.ones(len(genes), dtype=bool)
            if 'gene_n_cells' in z:
                coverage = np.asarray(z['gene_n_cells'][i])
                n = int(stars.n_cells.iloc[i])
                if coverage.shape != (len(genes),) or (coverage < 0).any() or (coverage > n).any():
                    raise ValueError('Invalid per-gene observed-cell support')
                expected_finite = coverage == n
            if (supported[i] and (not np.isfinite(values[expected_finite]).all() or not np.isnan(values[~expected_finite]).all())) or (not supported[i] and not np.isnan(values).all()):
                raise ValueError('Invalid missingness or finite-value convention')
    z.close()
    for sid, group in stars.groupby('sample_id'):
        n = int(samples.loc[samples.sample_id.eq(sid), 'retained_cells'].iloc[0])
        if int(group.n_cells.sum()) != n:
            raise ValueError('Star counts do not reconcile to sample cells')
    return report

def process(config, annotation=None):
    """Process a reviewed configuration; cleanup happens only after validation."""
    cfg = json.loads(Path(config).read_text(encoding='utf-8')) if not isinstance(config, dict) else config
    pipeline_checksum = digest(__file__)
    np.random.seed(0)
    organ_scope = cfg.get('scope') == 'reviewed_organ'
    if organ_scope:
        if not cfg.get('reviewed_tissue') or not cfg.get('allowed_tissues') or cfg.get('cleanup_download', True):
            raise ValueError('Organ scope requires reviewed tissues and explicit source retention')
    elif not cfg.get('reviewed_pbmc'):
        raise ValueError('PBMC scope must be reviewed first')
    if not cfg.get('sample_columns') or not cfg.get('sample_mapping_evidence'):
        raise ValueError('Specimen mapping must be reviewed first')
    reviewed_donor_column = donor_identity_column(cfg)
    if cfg.get('annotation_normalization') != 'aligned_reference_tp10k_log1p':
        raise ValueError('Configuration must explicitly select the verified SCimilarity normalization')
    out = ROOT / 'results' / ('organ_star_maps' if organ_scope else 'star_maps') / cfg['dataset_id']
    if (out / 'run.json').exists():
        saved = validate_outputs(out)
        if saved['configuration'] != cfg:
            raise ValueError('Changed configuration requires a new output location/version')
        print('Already validated:', cfg['dataset_id'], flush=True)
        return saved
    # A corrected/versioned processing contract may intentionally reuse the
    # exact immutable carrier from an earlier key.  Keep this as an identifier,
    # not an arbitrary path, so sources remain confined to data/downloads.
    source_dataset_id = cfg.get('source_dataset_id', cfg['dataset_id'])
    if (not isinstance(source_dataset_id, str) or not source_dataset_id
            or Path(source_dataset_id).name != source_dataset_id):
        raise ValueError('source_dataset_id must be a single dataset identifier')
    source = ROOT / 'data' / 'downloads' / (source_dataset_id + '.h5ad')
    archive = None
    asset_paths = []
    if cfg.get('input_format') == 'h5ad_tar':
        import tarfile
        import shutil
        archive = ROOT / 'data' / 'downloads' / (cfg['dataset_id'] + '.h5ad.tar.gz')
        download(cfg['url'], archive, cfg['expected_bytes'])
        source_manifest = json.loads(archive.with_suffix(archive.suffix + '.manifest.json').read_text())
        with tarfile.open(archive) as handle:
            members = [m for m in handle if m.isfile() and m.name.endswith('.h5ad')]
            if len(members) != 1:
                raise ValueError('Expected exactly one H5AD archive member')
            member = members[0]
            if footprint() + member.size >= LIMIT or shutil.disk_usage(ROOT).free < member.size:
                raise RuntimeError('H5AD archive expansion exceeds available storage')
            # Copy only bytes to the fixed project path; never extract archive paths.
            with source.open('wb') as target:
                shutil.copyfileobj(handle.extractfile(member), target)
            source_manifest['selected_h5ad_member'] = member.name
    elif cfg.get('input_format') == 'pinned_h5_assignments':
        from src.pinned_h5_assignments import prepare
        source_manifest, asset_paths = prepare(cfg, source)
    elif cfg.get('input_format') == 'pinned_cell_assignments':
        from src.pinned_cell_assignments import prepare
        source_manifest, asset_paths = prepare(cfg, source)
    elif cfg.get('input_format') == 'pinned_author_counts':
        from src.pinned_author_counts import prepare
        source_manifest, asset_paths = prepare(cfg, source)
    elif cfg.get('input_format') == 'cmo_assignments':
        from src.cmo_assignments import prepare
        source_manifest, asset_paths = prepare(cfg, source)
    elif cfg.get('input_format') == 'atlas_geo_assignments':
        from src.atlas_geo_assignments import prepare
        source_manifest, asset_paths = prepare(cfg, source)
    elif cfg.get('input_format') == 'paired_pool_metadata':
        from src.paired_pool_metadata import prepare
        source_manifest, asset_paths = prepare(cfg, source)
    elif cfg.get('input_format') == 'pooled_pbmc_181897':
        from src.pooled_pbmc_181897 import prepare
        source_manifest, asset_paths = prepare(cfg, source)
    elif cfg.get('input_format') == 'gzip_h5ad':
        from src.gzip_h5ad import prepare
        source_manifest, asset_paths = prepare(cfg, source)
    elif cfg.get('input_format') == 'paired_mtx_tables':
        from src.paired_mtx_tables import prepare
        source_manifest, asset_paths = prepare(cfg, source)
    elif cfg.get('input_format') == 'rds_seurat':
        from src.rds_seurat import prepare
        source_manifest, asset_paths = prepare(cfg, source)
    elif cfg.get('input_format') == 'ivanova_original':
        from src.ivanova_original import prepare
        source_manifest, asset_paths = prepare(cfg, source)
    elif cfg.get('input_format') == 'native_mtx_bundle':
        from src.native_mtx_bundle import prepare
        source_manifest, asset_paths = prepare(cfg, source)
    elif cfg.get('input_format') == 'ircolitis_original':
        from src.ircolitis_original import prepare
        source_manifest, asset_paths = prepare(cfg, source)
    elif cfg.get('input_format') == 'rag_pbmc_original':
        from src.rag_pbmc_original import prepare
        source_manifest, asset_paths = prepare(cfg, source)
    elif cfg.get('input_format') == 'lymphoma_original':
        from src.lymphoma_original import prepare
        source_manifest, asset_paths = prepare(cfg, source)
    elif cfg.get('input_format') == 'cillo_original':
        from src.cillo_original import prepare
        source_manifest, asset_paths = prepare(cfg, source)
    elif cfg.get('input_format') == 'reviewed_native_bundle':
        from src.reviewed_native_bundle import prepare
        source_manifest, asset_paths = prepare(cfg, source)
    elif cfg.get('input_format') == 'reviewed_public_native_bundle':
        from src.public_native_bundle import prepare
        source_manifest, asset_paths = prepare(cfg, source)
    elif cfg.get('input_format') == 'memory_b_original':
        from src.memory_b_original import prepare
        source_manifest, asset_paths = prepare(cfg, source)
    elif cfg.get('input_format') == 'sorted_memory_original':
        from src.sorted_memory_original import prepare
        source_manifest, asset_paths = prepare(cfg, source)
    elif cfg.get('input_format') == 'monocyte_plate_original':
        from src.monocyte_plate_original import prepare
        source_manifest, asset_paths = prepare(cfg, source)
    elif cfg.get('input_format') == 'hsct_pbmc_original':
        from src.hsct_pbmc_original import prepare
        source_manifest, asset_paths = prepare(cfg, source)
    elif cfg.get('input_format') == 'diabetes_initial':
        from src.diabetes_initial import prepare
        source_manifest, asset_paths = prepare(cfg, source)
    elif cfg.get('input_format') == 'lyme_vireo':
        from src.lyme_vireo import prepare
        source_manifest, asset_paths = prepare(cfg, source)
    elif cfg.get('input_format') == 'critical_illness_raw':
        from src.critical_illness_raw import prepare
        source_manifest, asset_paths = prepare(cfg, source)
    elif cfg.get('input_format') == 'partitioned_h5ad_counts':
        from src.partitioned_h5ad_counts import prepare
        source_manifest, asset_paths = prepare(cfg, source)
    elif cfg.get('input_format') == 'count_layer_h5ad':
        from src.count_layer_h5ad import prepare
        source_manifest, asset_paths = prepare(cfg, source)
    elif cfg.get('input_format') == 'typed_h5ad':
        from src.typed_h5ad import prepare
        source_manifest, asset_paths = prepare(cfg, source)
    elif cfg.get('input_format') == 'rdata_counts':
        from src.rdata_counts import prepare
        source_manifest, asset_paths = prepare(cfg, source)
    elif cfg.get('input_format') == 'hca_eqtl45':
        from src.hca_eqtl45 import prepare
        source_manifest, asset_paths = prepare(cfg, source)
    elif cfg.get('input_format') == 'influenza_hashing':
        from src.influenza_hashing import prepare
        source_manifest, asset_paths = prepare(cfg, source)
    elif cfg.get('input_format') == 'cll_hashing':
        from src.cll_hashing import prepare
        source_manifest, asset_paths = prepare(cfg, source)
    elif cfg.get('input_format') == 'bcell_hashing':
        from src.bcell_hashing import prepare
        source_manifest, asset_paths = prepare(cfg, source)
    elif cfg.get('input_format') == 'molecule_info':
        from src.molecule_info import prepare
        source_manifest, asset_paths = prepare(cfg, source)
    elif cfg.get('input_format') == 'reviewed_panel_h5':
        from src.reviewed_panel_h5 import prepare
        source_manifest, asset_paths = prepare(cfg, source)
    elif cfg.get('input_format') == 'prefix_mtx':
        from src.prefix_mtx import prepare
        source_manifest, asset_paths = prepare(cfg, source)
    elif cfg.get('input_format') == 'benchmark_archive_RNA':
        from src.benchmark_archive_RNA import prepare
        source_manifest, asset_paths = prepare(cfg, source)
    elif cfg.get('input_format') == 'parse_geo':
        from src.parse_geo import prepare
        source_manifest, asset_paths = prepare(cfg, source)
    elif cfg.get('input_format') == 'reviewed_series_mtx':
        from src.reviewed_series_mtx import prepare
        source_manifest, asset_paths = prepare(cfg, source)
    elif cfg.get('input_format') == 'geo_mtx_table':
        from src.geo_mtx_table import prepare
        source_manifest, asset_paths = prepare(cfg, source)
    elif cfg.get('input_format') == 'geo_10x_tar':
        from src.geo_inputs import convert_geo_tar
        archive = ROOT / 'data' / 'downloads' / (cfg['dataset_id'] + '_RAW.tar')
        download(cfg['url'], archive, cfg['expected_bytes'])
        source_manifest = json.loads(archive.with_suffix('.tar.manifest.json').read_text())
        convert_geo_tar(archive, source, cfg)
    elif cfg.get('input_format') == 'geo_10x_assets':
        from src.geo_assets import prepare_assets
        from src.geo_inputs import convert_geo_tar
        archive, source_manifest, asset_paths = prepare_assets(cfg)
        convert_geo_tar(archive, source, cfg)
    else:
        download(cfg['url'], source, cfg['expected_bytes'])
        source_manifest = json.loads(source.with_suffix('.h5ad.manifest.json').read_text())
    a = ad.read_h5ad(source, backed='r')
    try:
        from src.source_metadata import apply_metadata_mapping
        obs = apply_metadata_mapping(a.obs, cfg.get('metadata_mapping', {}))
        for column in dict.fromkeys(cfg['sample_columns'] + [reviewed_donor_column]):
            if column not in obs or obs[column].isna().any() or obs[column].astype(str).str.lower().isin(['unknown', 'nan', 'na', '']).any():
                raise ValueError('Missing or unknown specimen identity: ' + column)
        if cfg['raw_location'] == 'raw.X' and a.raw is not None:
            count_matrix = a.raw.X; var = a.raw.var.copy()
        elif cfg['raw_location'] == 'X':
            count_matrix = a.X; var = a.var.copy()
        else:
            raise ValueError('An explicitly verified raw count location is required')
        availability = None
        source_libraries = None
        if 'source_gene_availability' in a.uns and not cfg.get('use_source_gene_availability', False):
            raise ValueError('Source gene availability requires explicit configuration')
        if cfg.get('use_source_gene_availability', False):
            from src.gene_availability import validate_availability
            if not np.array_equal(np.asarray(a.uns['source_gene_ids']).astype(str), var.index.astype(str).to_numpy()):
                raise ValueError('Source gene availability does not match gene order')
            source_libraries, availability = validate_availability(obs['source_library_index'].to_numpy(), a.uns['source_gene_availability'], len(var))
        symbols = annotation_symbols(var, cfg)
        output_source_indices = np.arange(len(var), dtype=np.int64)
        output_var = var.copy()
        output_panel_path = cfg.get('star_output_gene_panel')
        if output_panel_path is not None:
            reason = cfg.get('star_output_gene_panel_reason')
            if not isinstance(reason, str) or not reason.strip():
                raise ValueError('A bounded star-output panel requires an explicit reason')
            panel_path = ROOT / str(output_panel_path)
            panel = [line.strip() for line in panel_path.read_text(encoding='utf-8').splitlines() if line.strip()]
            from src.star_storage import output_gene_indices
            output_source_indices = output_gene_indices(symbols, panel)
            output_var = var.iloc[output_source_indices].copy()
            output_var['source_gene_index'] = output_source_indices
        from src.source_selection import primary_data_mask
        selected = primary_data_mask(obs, cfg.get('component_primary_override'))
        allowed_tissues = cfg['allowed_tissues'] if organ_scope else ['blood']
        selected &= obs['tissue'].astype(str).isin(allowed_tissues).to_numpy()
        from src.source_selection import source_exclusion_mask
        source_excluded = source_exclusion_mask(obs, cfg.get('exclude_source_values', {}))
        selected &= ~source_excluded
        obs['collection_source_excluded'] = source_excluded
        obs['sample_id'] = cfg['study_id'] + ':' + obs[cfg['sample_columns']].astype(str).agg('|'.join, axis=1)
        obs['donor_uid'] = cfg['study_id'] + ':' + obs[reviewed_donor_column].astype(str)
        # Prevent a composite specimen key from accidentally grouping different donors.
        if obs.groupby('sample_id', observed=True).donor_uid.nunique().max() > 1:
            raise ValueError('A sample key contains several donors')
        from src.specimen_views import attach_physical_identity
        obs=attach_physical_identity(obs,cfg)
        qc = np.zeros(a.n_obs, dtype=bool)
        libraries = np.zeros(a.n_obs)
        detected = np.zeros(a.n_obs, dtype=int)
        mito = np.zeros(a.n_obs)
        mito_genes = np.char.startswith(symbols.astype(str), 'MT-')
        mito_unavailable = cfg.get('qc_mito_unavailable', False)
        if not isinstance(mito_unavailable, bool):
            raise ValueError('qc_mito_unavailable must be an explicit boolean')
        if mito_unavailable:
            reason = cfg.get('qc_mito_unavailable_reason')
            if not isinstance(reason, str) or not reason.strip():
                raise ValueError('Unavailable mitochondrial QC requires a source-specific reason')
            if mito_genes.any() or cfg.get('qc_mito_percent_column'):
                raise ValueError('Unavailable mitochondrial QC conflicts with native mitochondrial genes or source percentages')
        block_size = cfg.get('block_size', 1024)
        for start in range(0, a.n_obs, block_size):
            stop = min(start + block_size, a.n_obs)
            x = sparse.csr_matrix(count_matrix[start:stop], dtype=np.float64)
            if availability is not None:
                from src.gene_availability import require_observed_counts
                require_observed_counts(x, source_libraries[start:stop], availability)
            from src.count_semantics import validate_count_values
            validate_count_values(x.data,cfg)
            x.eliminate_zeros()
            lib = np.asarray(x.sum(axis=1)).ravel()
            genes = np.diff(x.indptr)
            fraction = np.divide(np.asarray(x[:, mito_genes].sum(axis=1)).ravel(), lib, out=np.zeros_like(lib), where=lib > 0)
            if cfg.get('qc_mito_percent_column'):
                source_percent = pd.to_numeric(obs[cfg['qc_mito_percent_column']].iloc[start:stop], errors='raise').to_numpy()
                if not np.isfinite(source_percent).all() or (source_percent < 0).any() or (source_percent > 100).any():
                    raise ValueError('Invalid source mitochondrial percentage')
                fraction = source_percent / 100
            if mito_unavailable:
                fraction = np.full_like(lib, np.nan)
            libraries[start:stop] = lib; detected[start:stop] = genes; mito[start:stop] = fraction
            mito_gate = np.ones(len(lib), dtype=bool) if mito_unavailable else fraction <= cfg['qc']['max_mito_fraction']
            qc[start:stop] = (lib >= cfg['qc']['min_counts']) & (genes >= cfg['qc']['min_genes']) & mito_gate & (lib > 0)
        retained = selected & qc
        if not retained.any():
            raise ValueError('No cells pass the reviewed selection and QC')
        if annotation is None:
            import torch
            from scimilarity import CellAnnotation
            torch.set_num_threads(8)
            torch.manual_seed(0)
            print('Loading SCimilarity reference', flush=True)
            annotation = CellAnnotation(str(ROOT / 'models' / 'model_v1.1'), use_gpu=torch.cuda.is_available())
        mapper = gene_alignment(symbols, annotation.gene_order)
        if int((np.asarray(mapper.sum(axis=0)).ravel() > 0).sum()) < 5000:
            raise ValueError('Fewer than 5,000 reference genes overlap; SCimilarity default compatibility threshold is not met')
        labels = np.full(a.n_obs, '', dtype=object)
        vote = np.full(a.n_obs, np.nan); distance = vote.copy()
        annotation_libraries = np.full(a.n_obs, np.nan)
        kept_indices = np.flatnonzero(retained)
        for start in range(0, len(kept_indices), block_size):
            idx = kept_indices[start:start + block_size]
            x = sparse.csr_matrix(count_matrix[idx], dtype=np.float64)
            # Annotation follows the fixed reference's documented preprocessing;
            # expression-star normalization separately uses all source RNA genes.
            aligned, reference_library = annotation_input(x, mapper)
            annotation_libraries[idx] = reference_library
            embedding = annotation.get_embeddings(aligned.astype(np.float32), buffer_size=block_size)
            predictions, _, _, stats = annotation.get_predictions_knn(embedding, k=50, ef=100, weighting=False, disable_progress=True)
            labels[idx] = predictions.to_numpy(); vote[idx] = stats.vsAll.to_numpy(); distance[idx] = stats.min_dist.to_numpy()
            print(f"Annotated {min(start + block_size, len(kept_indices))}/{len(kept_indices)} cells", flush=True)
        out.mkdir(parents=True, exist_ok=True)
        # Unwritten NaN-filled rows do not allocate chunks. Budget actual
        # supported rows, using labels already computed, without a dense fiction.
        n_samples = obs.loc[selected, 'sample_id'].nunique()
        from src.star_storage import export_budget
        output_keys = output_statistic_keys(cfg)
        storage_plan=export_budget(labels[retained],obs.loc[retained,'sample_id'],cfg['min_cells_per_star'],len(output_var),len(annotation.classes),
            obs_bytes=int(obs.memory_usage(deep=True).sum()),availability_bytes=int(availability.nbytes) if availability is not None else 0,view_count=n_samples,
            array_count=len(output_keys))
        expected_output=storage_plan['expected_bytes']
        if footprint() + expected_output >= LIMIT:
            raise RuntimeError('Star export would exceed project disk limit')
        rows = []
        expression = h5py.File(out / 'expression.h5', 'w')
        arrays = {key: expression.create_dataset(key, shape=(n_samples * len(annotation.classes), len(output_var)),
            dtype='float64', chunks=(1, len(output_var)), compression='gzip', fillvalue=np.nan) for key in output_keys}
        gene_support = None
        if availability is not None:
            gene_support = expression.create_dataset('gene_n_cells', shape=(n_samples*len(annotation.classes),len(output_var)), dtype='int64', chunks=(1,len(output_var)), compression='gzip', fillvalue=0)
            expression.create_dataset('source_gene_availability', data=availability[:, output_source_indices], compression='gzip')
            expression.create_dataset('source_library_names', data=np.asarray(a.uns['source_library_names'],dtype=object), dtype=h5py.string_dtype('utf-8'))
        sample_rows = []
        confusion = []
        all_labels = sorted(annotation.classes)
        mapping = pd.read_csv(ROOT / 'references' / 'vocabulary' / 'label_mapping.tsv', sep='\t')
        for sid in sorted(obs.loc[selected, 'sample_id'].unique()):
            in_sample = obs.sample_id.eq(sid).to_numpy()
            sample_idx = np.flatnonzero(in_sample & retained)
            # Read a specimen once; repeated sparse HDF5 gathers per label can be costly.
            sample_counts = sparse.csr_matrix(count_matrix[sample_idx], dtype=np.float64) if len(sample_idx) else sparse.csr_matrix((0, len(var)))
            original = obs.loc[in_sample & selected].copy()
            meta = dict(sample_id=sid, donor_uid=original.donor_uid.iloc[0], source_donor_id=str(original[reviewed_donor_column].iloc[0]),
                source_cells=len(original), retained_cells=len(sample_idx), excluded_cells=int(len(original) - len(sample_idx)))
            # Store constant source fields; never choose the first value of conflicting metadata.
            for column in original.columns:
                unique = original[column].dropna().astype(str).unique()
                if len(unique) == 1 and column not in meta:
                    meta[column] = unique[0]
            meta['sample_mapping_evidence'] = cfg['sample_mapping_evidence']
            sample_rows.append(meta)
            for label in all_labels:
                idx = sample_idx[labels[sample_idx] == label]
                local_idx = np.flatnonzero(labels[sample_idx] == label)
                n = len(idx); supported = n >= cfg['min_cells_per_star']
                row = dict(sample_id=sid, label=label, n_cells=n,
                    fraction_of_retained=n / len(sample_idx) if len(sample_idx) else np.nan,
                    status='supported' if supported else ('absent' if n == 0 else 'insufficient_support'),
                    mean_library_size=float(libraries[idx].mean()) if n else np.nan,
                    total_raw_counts=float(libraries[idx].sum()),
                    mean_vote_fraction=float(vote[idx].mean()) if n else np.nan,
                    mean_min_neighbor_distance=float(distance[idx].mean()) if n else np.nan)
                output_index = len(rows)
                rows.append(row)
                if gene_support is not None and n:
                    from src.gene_availability import observed_gene_support
                    gene_support[output_index] = observed_gene_support(source_libraries[idx],availability,len(var))[output_source_indices]
                if supported:
                    x = sample_counts[local_idx]
                    for rep, matrix in zip(REPRESENTATIONS, representations(x)):
                        mean, variance, sem, total = summarize(matrix)
                        if availability is not None:
                            from src.gene_availability import mask_incomplete_summary
                            (mean,variance,sem,total), _ = mask_incomplete_summary((mean,variance,sem,total),source_libraries[idx],availability)
                        mean = mean[output_source_indices]
                        variance = variance[output_source_indices]
                        sem = sem[output_source_indices]
                        total = total[output_source_indices]
                        for name, values in [('mean', mean), ('variance', variance), ('sem', sem)]:
                            key = rep + '_' + name
                            if key in arrays:
                                arrays[key][output_index] = values
                        if rep == 'raw' and 'raw_sum' in arrays:
                            arrays['raw_sum'][output_index] = total
                else:
                    for key in arrays:
                        if key == 'raw_sum' and n:
                            values = np.asarray(sample_counts[local_idx].sum(axis=0)).ravel()[output_source_indices]
                            if availability is not None:
                                values[observed_gene_support(source_libraries[idx],availability,len(var))[output_source_indices] != n] = np.nan
                            arrays[key][output_index] = values
            if len(sample_idx):
                comparison = pd.DataFrame({'scimilarity_label': labels[sample_idx], 'author_cell_type': obs.iloc[sample_idx]['author_cell_type'].astype(str).to_numpy() if 'author_cell_type' in obs else 'not_reported'})
                table = comparison.groupby(['scimilarity_label', 'author_cell_type']).size().reset_index(name='n_cells')
                table.insert(0, 'sample_id', sid); confusion.append(table)
        pd.DataFrame(rows).to_csv(out / 'stars.csv', index=False)
        pd.DataFrame(sample_rows).to_csv(out / 'samples.csv', index=False)
        pd.DataFrame(sample_rows)[['donor_uid', 'source_donor_id']].drop_duplicates().to_csv(out / 'donors.csv', index=False)
        gene_table = output_var.copy()
        if 'source_gene_id' in gene_table.columns:
            # Some source-native H5ADs preserve the exact feature ID both as
            # the var index and as an explicit provenance column.  Permit the
            # redundant representation only when it is losslessly identical;
            # otherwise resetting the index would either fail or silently
            # choose between conflicting IDs.
            if not np.array_equal(
                gene_table['source_gene_id'].astype(str).to_numpy(),
                gene_table.index.astype(str).to_numpy(),
            ):
                raise ValueError('var source_gene_id column conflicts with the source gene index')
            gene_table = gene_table.drop(columns=['source_gene_id'])
        gene_table.rename_axis('source_gene_id').reset_index().to_csv(out / 'genes.csv', index=False)
        expression.close()
        # A small per-cell metadata audit enables tracing labels and QC without retaining expression matrices.
        audit = obs.copy(); audit['selected'] = selected; audit['qc_pass'] = qc; audit['retained'] = retained
        audit['scimilarity_label'] = labels; audit['library_size'] = libraries; audit['detected_genes'] = detected
        audit['mito_fraction'] = mito; audit['vote_fraction'] = vote; audit['min_neighbor_distance'] = distance
        if mito_unavailable:
            audit['mitochondrial_qc_available'] = False
        audit['annotation_reference_library_size'] = annotation_libraries
        audit.to_csv(out / 'cell_metadata.csv.gz', compression='gzip')
        if confusion:
            pd.concat(confusion, ignore_index=True).to_csv(out / 'annotation_comparison.csv', index=False)
        mapping.to_csv(out / 'annotation_label_mapping.tsv', sep='\t', index=False)
        (out / 'source_manifest.json').write_text(json.dumps(source_manifest, indent=2), encoding='utf-8')
        (out / 'dataset_source.json').write_text(json.dumps(cfg['source_metadata'], indent=2), encoding='utf-8')
        if 'source_assignment_audit' in a.uns:
            audit = {library: {key: int(value) for key, value in values.items()}
                     for library, values in a.uns['source_assignment_audit'].items()}
            (out / 'source_assignment_audit.json').write_text(json.dumps(audit, indent=2), encoding='utf-8')
        checksums = {p.name: digest(p) for p in out.iterdir() if p.is_file() and p.name != 'run.json'}
        report = dict(status='validated', completed_at=datetime.now(timezone.utc).isoformat(), configuration=cfg,
            cells_retained=int(retained.sum()), biological_samples=n_samples, donors=len(set(obs.loc[selected, 'donor_uid'])),
            annotation_gene_overlap=int((np.asarray(mapper.sum(axis=0)).ravel() > 0).sum()), reference_genes=len(annotation.gene_order),
            versions={p: importlib.metadata.version(p) for p in ('scimilarity', 'torch', 'anndata', 'numpy', 'scipy', 'pandas')},
            model_release='v1.1', model_provenance=json.loads((ROOT / 'environment' / 'model_manifest.json').read_text()),
            encoder_sha256=digest(ROOT / 'models' / 'model_v1.1' / 'encoder.ckpt'),
            gene_order_sha256=digest(ROOT / 'models' / 'model_v1.1' / 'gene_order.tsv'),
            pipeline_sha256=pipeline_checksum,
            auxiliary_code_sha256={name: digest(ROOT / name) for name in ['src/geo_inputs.py', 'src/geo_assets.py', 'src/source_metadata.py', 'src/source_selection.py', 'scripts/download_pbmc_asset.py'] + (['src/influenza_hashing.py','scripts/audit_influenza_hashing.py','src/geo_mtx_table.py'] if cfg.get('input_format') == 'influenza_hashing' else []) + (['src/molecule_info.py'] if cfg.get('input_format') == 'molecule_info' else []) + (['src/prefix_mtx.py','src/geo_mtx_table.py'] if cfg.get('input_format') == 'prefix_mtx' else []) + (['src/geo_mtx_table.py'] if cfg.get('input_format') == 'geo_mtx_table' else []) + (['src/hca_eqtl45.py'] if cfg.get('input_format') == 'hca_eqtl45' else []) + (['src/rdata_counts.py'] if cfg.get('input_format') == 'rdata_counts' else []) + (['src/typed_h5ad.py'] if cfg.get('input_format') == 'typed_h5ad' else []) + (['src/rds_seurat.py', 'src/rdata_counts.py'] + (['src/rds_streaming.py'] if 'rds_compression_layers' in cfg else []) if cfg.get('input_format') == 'rds_seurat' else []) + (['src/gzip_h5ad.py'] if cfg.get('input_format') == 'gzip_h5ad' else []) + (['src/pinned_cell_assignments.py'] if cfg.get('input_format') == 'pinned_cell_assignments' else []) + (['src/cmo_assignments.py'] if cfg.get('input_format') == 'cmo_assignments' else []) + (['src/atlas_geo_assignments.py'] if cfg.get('input_format') == 'atlas_geo_assignments' else []) + (['src/paired_pool_metadata.py'] if cfg.get('input_format') == 'paired_pool_metadata' else [])},
            output_sha256=checksums,
            uncertainty=(
                ('Per-gene unbiased cell variance and SEM=sqrt(variance/n) are persisted; '
                 if any(key.endswith('_variance') for key in output_keys)
                 else 'Per-gene SEM is persisted and unbiased cell variance is exactly recoverable as SEM squared times n_cells; ')
                + 'conditional on QC and annotation, not donor replication or annotation uncertainty.'
                + (' Gene statistics require all contributing cells to have a source-available feature; incomplete genes are NaN with observed-cell support retained.' if availability is not None else '')
            ),
            annotation_confidence='Uncalibrated reference vote fraction and nearest-neighbor distance; no confidence-based exclusion or PBMC label restriction.',
            source_preconversion_audit=json.loads(a.uns.get('source_preconversion_audit_json', '{}')),
            processing_generation=2,
            normalization=dict(annotation='Raw counts aligned and duplicate symbols summed in reference gene order, then TP10K normalized over aligned genes and natural-log log1p.',
                stars='Raw mean; mean of per-cell TP10K over all source RNA genes; mean of cell-level natural-log log1p(TP10K).'),
            seed=0, annotation_settings=dict(k=50, ef=100, weighting=False, label_restriction=None),
            limitations=['No new doublet detection; source-filtered cells plus provisional QC.', 'No cross-study batch correction.', 'No discovery or donor splits assigned yet.', 'Sample and donor counts are study-scoped; cross-source deduplication is still required.'])
        if mito_unavailable:
            report['mitochondrial_qc'] = dict(evaluated=False, available=False,
                native_mitochondrial_symbol_features=0, reason=cfg['qc_mito_unavailable_reason'],
                policy='Mitochondrial fractions are NaN; only the minimum RNA count, detected-gene and positive-library QC gates are evaluated.')
            report['limitations'].append('Mitochondrial QC is unavailable for this source feature axis; a zero fraction is not measured or imputed.')
        if cfg.get('physical_sample_columns') is not None:
            report['observation_views']=n_samples
            report['biological_samples']=int(obs.loc[selected,'physical_specimen_uid'].nunique())
            report['physical_sample_mapping_evidence']=cfg['physical_sample_mapping_evidence']
            report['auxiliary_code_sha256']['src/specimen_views.py']=digest(ROOT/'src/specimen_views.py')
            report['limitations'].append('Centroids are conditional observation views; do not pool sorting views into a physical-specimen map or treat them as independent specimens.')
        if cfg.get('input_format') == 'count_layer_h5ad':
            report['auxiliary_code_sha256']['src/count_layer_h5ad.py']=digest(ROOT/'src/count_layer_h5ad.py')
        if cfg.get('input_format') == 'partitioned_h5ad_counts':
            report['auxiliary_code_sha256']['src/partitioned_h5ad_counts.py']=digest(ROOT/'src/partitioned_h5ad_counts.py')
        if cfg.get('input_format') == 'pinned_h5_assignments':
            report['auxiliary_code_sha256']['src/pinned_h5_assignments.py']=digest(ROOT/'src/pinned_h5_assignments.py')
        if cfg.get('input_format') == 'paired_mtx_tables':
            for name in ['src/paired_mtx_tables.py', 'src/geo_mtx_table.py']:
                report['auxiliary_code_sha256'][name] = digest(ROOT / name)
        if any('host_gene_selection' in spec for spec in cfg.get('geo_samples', [])):
            name = 'src/human_ensembl_rna.py'
            report['auxiliary_code_sha256'][name] = digest(ROOT / name)
        if cfg.get('input_format') == 'pinned_author_counts':
            report['auxiliary_code_sha256']['src/pinned_author_counts.py']=digest(ROOT/'src/pinned_author_counts.py')
        report['star_storage_budget']=storage_plan
        report['persisted_star_statistics']=list(output_keys)
        if tuple(output_keys) != DEFAULT_OUTPUT_KEYS:
            report['limitations'].append(
                'Persisted star statistics are an explicit bounded subset; the complete retained raw-count carrier remains the source for replay.'
            )
        if cfg.get('input_format') == 'benchmark_archive_RNA':
            report['auxiliary_code_sha256']['src/benchmark_archive_RNA.py']=digest(ROOT/'src/benchmark_archive_RNA.py')
        if cfg.get('input_format') == 'pooled_pbmc_181897':
            report['auxiliary_code_sha256']['src/pooled_pbmc_181897.py']=digest(ROOT/'src/pooled_pbmc_181897.py')
        report['auxiliary_code_sha256']['src/star_storage.py']=digest(ROOT/'src/star_storage.py')
        if any(spec.get('rna_table_layout') == 'three_gene_identifiers' for spec in cfg.get('geo_samples', [])):
            report['auxiliary_code_sha256']['src/annotated_gene_table.py'] = digest(ROOT/'src/annotated_gene_table.py')
        if any(spec.get('stream_rna_table') is True for spec in cfg.get('geo_samples', [])):
            report['auxiliary_code_sha256']['src/streamed_gene_table.py'] = digest(ROOT/'src/streamed_gene_table.py')
        if any('gene_feature_selection' in v for v in cfg.get('geo_samples',[])):
            report['auxiliary_code_sha256']['src/gene_te_h5.py']=digest(ROOT/'src/gene_te_h5.py')
        if availability is not None:
            report['auxiliary_code_sha256']['src/gene_availability.py'] = digest(ROOT/'src/gene_availability.py')
            report['source_gene_availability'] = 'Per-library observed genes and per-star gene_n_cells retained in expression.h5. An all-cell gene statistic is undefined when observed support differs from n_cells.'
            report['normalization']['stars'] = 'Raw mean; per-cell TP10K over observed source RNA genes before averaging; cell-level log1p(TP10K) before averaging. Incomplete genes remain NaN; zero union placeholders never count as observed genes.'
        if cfg.get('input_format') == 'parse_geo':
            for name in ['src/parse_geo.py', 'scripts/audit_parse_transplant_coordinates.py', 'scripts/prepare_parse_transplant_295359.py']:
                report['auxiliary_code_sha256'][name] = digest(ROOT/name)
        if cfg.get('input_format') == 'cll_hashing':
            for name in ['src/cll_hashing.py', 'src/hashsolo_assignments.py', 'scripts/audit_cll_hashsolo_283397.py', 'scripts/review_cll_hashing_283397.py']:
                report['auxiliary_code_sha256'][name] = digest(ROOT/name)
            report['limitations'][0] = 'Hash-based doublet/negative/low-confidence exclusions; no RNA-based or same-specimen doublet detection.'
        if cfg.get('input_format') == 'bcell_hashing':
            for name in ['src/bcell_hashing.py','src/hashsolo_assignments.py','scripts/audit_bcell_hashsolo_207475.py','scripts/audit_bcell_pools_207475.py']:
                report['auxiliary_code_sha256'][name] = digest(ROOT/name)
            report['limitations'][0] = 'Hash-based doublet/negative/low-confidence exclusions and observed RNA/tag intersection; no RNA-based or same-specimen doublet detection.'
        if cfg.get('input_format') == 'diabetes_initial':
            for name in ['src/diabetes_initial.py','src/rds_seurat.py','src/rdata_counts.py']:
                report['auxiliary_code_sha256'][name] = digest(ROOT/name)
        if cfg.get('input_format') == 'lyme_vireo':
            report['auxiliary_code_sha256']['src/lyme_vireo.py'] = digest(ROOT/'src/lyme_vireo.py')
        if cfg.get('input_format') == 'critical_illness_raw':
            for name in ['src/critical_illness_raw.py','src/cellranger_called_cells.py']:
                report['auxiliary_code_sha256'][name] = digest(ROOT/name)
        if cfg.get('input_format') == 'ivanova_original':
            for name in ['src/ivanova_original.py', 'src/rdata_counts.py']:
                report['auxiliary_code_sha256'][name] = digest(ROOT/name)
            report['limitations'].append('Exact proved original RNA copies excluded by a pinned map. Unresolved LB material labels are exported separately and excluded from acquired specimen counts.')
        if cfg.get('input_format') == 'native_mtx_bundle':
            report['auxiliary_code_sha256']['src/native_mtx_bundle.py'] = digest(ROOT/'src/native_mtx_bundle.py')
            report['limitations'].append('Physical specimen UIDs are conservative patient/day accounting groups. Exact shared tubes/aliquots are unverified; captures remain separate centroids and are not pooled. Existing atlas coordinates are not silently remapped.')
        if cfg.get('input_format') == 'ircolitis_original':
            for name in ['src/ircolitis_original.py', 'scripts/audit_PBMC_206298_original_counts.py']:
                report['auxiliary_code_sha256'][name] = digest(ROOT/name)
            report['limitations'].append('Initial-only public source objects; two GEO followup records absent. Source day is not a therapy clock. Capture selection and original GEO tissue discrepancies remain explicit; exact repeated source cell counted once with both author cluster labels retained.')
        if cfg.get('input_format') == 'rag_pbmc_original':
            for name in ['src/rag_pbmc_original.py', 'scripts/audit_PBMC_274284_original_counts.py', 'scripts/review_PBMC_274284_assignments.py']:
                report['auxiliary_code_sha256'][name] = digest(ROOT/name)
            report['limitations'].append('Author SNG/log-likelihood selection is not calibrated assignment probability. Raw whitelist rows are not cells; bulk-record ages are not assumed scRNA draw ages. Nativehg19gene symbols retained without modernization; sorted and unsorted capture views remain dependent.')
        if cfg.get('input_format') == 'lymphoma_original':
            for name in ['src/lymphoma_original.py', 'src/rag_pbmc_original.py', 'scripts/prepare_PBMC_180578.py', 'scripts/prepare_PBMC_207192.py', 'scripts/audit_PBMC_207192_original_post_counts.py']:
                report['auxiliary_code_sha256'][name] = digest(ROOT/name)
            report['limitations'].append(cfg['source_metadata']['source_limitations'])
        if cfg.get('input_format') == 'cillo_original':
            for name in ['src/cillo_original.py', 'src/rag_pbmc_original.py', 'scripts/prepare_PBMC_180578.py', 'scripts/stage_PBMC_180578_native_counts.py', 'scripts/review_PBMC_180578_specimens.py', 'scripts/audit_PBMC_180578_original_count_identity.py']:
                report['auxiliary_code_sha256'][name] = digest(ROOT/name)
            report['limitations'].append(cfg['source_metadata']['source_limitations'])
        if cfg.get('input_format') == 'reviewed_native_bundle':
            for name in ['src/reviewed_native_bundle.py', 'src/rag_pbmc_original.py', 'scripts/audit_PBMC_' + cfg['study_id'][3:] + '_original_counts.py']:
                report['auxiliary_code_sha256'][name] = digest(ROOT/name)
            report['limitations'].append(cfg['source_metadata']['source_limitations'])
        if cfg.get('input_format') == 'reviewed_public_native_bundle':
            for name in ['src/public_native_bundle.py', 'src/rag_pbmc_original.py', cfg['native_replay_code']['path']] + cfg['native_source_code_paths']:
                report['auxiliary_code_sha256'][name] = digest(ROOT/name)
            report['limitations'].append(cfg['source_metadata']['source_limitations'])
        if cfg.get('input_format') == 'memory_b_original':
            for name in ['src/memory_b_original.py', 'src/rag_pbmc_original.py', 'scripts/audit_PBMC_157447_original_counts.py']:
                report['auxiliary_code_sha256'][name] = digest(ROOT/name)
            report['limitations'].append('All825provided source RNAcells before ownQC; author smaller plate/Ig-gene analysis selection is not imposed. Literal plate flag matches132source cells and does not reproduce stated180cell analysis. FU/LA labels remain opaque;20dependent views across12conservative donor/material groups. Original batch concerns remain; not whole-PBMC composition.')
        if cfg.get('input_format') == 'sorted_memory_original':
            for name in ['src/sorted_memory_original.py', 'src/rag_pbmc_original.py', 'scripts/audit_PBMC_106540_original_counts.py']:
                report['auxiliary_code_sha256'][name] = digest(ROOT/name)
            report['limitations'].append('Sorted PBMC-derived CD4memory subsets are dependent views across16conservative participant/material groups, not established distinct blood draws. Infection history is not current disease. Exact GENCODE19versions/symbols retained; HTSeq read counts are not relabelled UMIs. No whole-PBMC composition claim.')
        if cfg.get('input_format') == 'monocyte_plate_original':
            for name in ['src/monocyte_plate_original.py', 'src/rag_pbmc_original.py', 'scripts/audit_PBMC_184241_original_counts.py']:
                report['auxiliary_code_sha256'][name] = digest(ROOT/name)
            report['limitations'].append('Isolated CD14monocyte source;18ex vivo views across9donor/visit groups and3donors. Missing mitochondrial coordinates make mitochondrial QC unavailable. Literal workbook-ID join is independent of workbook row order; no whole-PBMC composition claim.')
        if cfg.get('input_format') == 'hsct_pbmc_original':
            for name in ['src/hsct_pbmc_original.py', 'src/rag_pbmc_original.py', 'scripts/audit_PBMC_266519_original_counts.py']:
                report['auxiliary_code_sha256'][name] = digest(ROOT/name)
            report['limitations'].append('CD3/CD19-enriched RNA subset further conditioned on explicit author identity annotation availability.47148source RNAcells remain unassigned; no identity inferred from numeric barcode suffix. Repeated visits are not independent donors; no whole-PBMC or unconditional B/T composition claim.')
        if cfg.get('input_format') == 'reviewed_panel_h5':
            report['auxiliary_code_sha256']['src/reviewed_panel_h5.py'] = digest(ROOT/'src/reviewed_panel_h5.py')
        (out / 'run.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
        validate_outputs(out)
    finally:
        a.file.close()
    # Delete only the exact temporary input after all saved output checks pass.
    if cfg.get('cleanup_download', True):
        resolved = source.resolve()
        if not resolved.is_relative_to((ROOT / 'data' / 'downloads').resolve()):
            raise ValueError('Unsafe cleanup path')
        resolved.unlink()
        if archive is not None:
            resolved_archive = archive.resolve()
            if not resolved_archive.is_relative_to((ROOT / 'data' / 'downloads').resolve()):
                raise ValueError('Unsafe archive cleanup path')
            resolved_archive.unlink()
        for asset in asset_paths:
            resolved_asset = asset.resolve()
            if not resolved_asset.is_relative_to((ROOT / 'data' / 'downloads').resolve()):
                raise ValueError('Unsafe source asset cleanup path')
            resolved_asset.unlink()
    print(json.dumps({k: report[k] for k in ('biological_samples', 'donors', 'cells_retained', 'annotation_gene_overlap')}, indent=2), flush=True)
    return report
