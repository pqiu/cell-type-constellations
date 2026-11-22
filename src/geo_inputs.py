"""Load reviewed per-specimen 10x matrices from a GEO tar without extracting paths."""
import gzip
import csv
import io
import json
import subprocess
import tempfile
import uuid
import tarfile
import urllib.request
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import io as scipy_io
from scipy import sparse

def soft_samples(text):
    result = {}; current = None
    for line in text.splitlines():
        if line.startswith('^SAMPLE = '):
            current = line.split(' = ', 1)[1]; result[current] = {}
        elif line.startswith('^'):
            current = None
        elif current and line.startswith('!Sample_') and ' = ' in line:
            key, value = line.split(' = ', 1)
            result[current].setdefault(key.removeprefix('!Sample_'), []).append(value)
    return result


def read_gene_table(raw, separator='\t'):
    """Read gene-by-cell GEX TSV, checking identity and count semantics."""
    text = raw.decode('utf-8')
    if separator not in ('\t', ',', ';', ' '):
        raise ValueError('Unsupported reviewed gene-table separator')
    lines = csv.reader(io.StringIO(text), delimiter=separator, skipinitialspace=separator == ' ')
    header = next(lines)
    first = next(lines)
    if len(first) == len(header) + 1:
        barcodes = header
    elif len(first) == len(header):
        barcodes = header[1:]
    else:
        raise ValueError('Unexpected GEX table header width')
    if not barcodes or len(set(barcodes)) != len(barcodes) or any(not b for b in barcodes):
        raise ValueError('GEX table cell identifiers must be unique and nonempty')
    blocks = []; names = []
    for frame in pd.read_csv(io.StringIO(text), sep=separator, index_col=0, chunksize=256,
                             keep_default_na=False, skipinitialspace=separator == ' '):
        if list(frame.columns) != barcodes:
            raise ValueError('GEX table cell order changed')
        x = frame.to_numpy(dtype=np.float64)
        if not np.isfinite(x).all() or (x < 0).any() or not np.equal(x, np.floor(x)).all():
            raise ValueError('GEX table does not contain nonnegative integer counts')
        names.extend(frame.index.astype(str))
        blocks.append(sparse.csr_matrix(x))
    if len(set(names)) != len(names):
        raise ValueError('Duplicate GEX gene identifiers')
    return sparse.vstack(blocks).T.tocsr(), pd.DataFrame({0: names, 1: names}), pd.Series(barcodes)

def restore_missing_dimensions(raw, n_genes, n_cells):
    """Explicit opt-in for a reviewed coordinate file lacking its dimension line."""
    lines = raw.splitlines()
    if not lines or lines[0].strip() != b'%%MatrixMarket matrix coordinate integer general':
        raise ValueError('Expected integer coordinate Matrix Market header')
    entries = []
    for line in lines[1:]:
        if not line.strip() or line.startswith(b'%'):
            continue
        values = line.split()
        if len(values) != 3:
            raise ValueError('Malformed coordinate entry')
        row, col, count = map(int, values)
        if not (1 <= row <= n_genes and 1 <= col <= n_cells and count >= 0):
            raise ValueError('Coordinate outside reviewed gene/barcode dimensions')
        entries.append(line)
    return (lines[0] + b'\n' + f'{n_genes} {n_cells} {len(entries)}\n'.encode()
            + b'\n'.join(entries) + b'\n')


def align_pool_metadata(raw, barcodes, assignment):
    """Strict source barcode join; never infer a donor for unassigned cells."""
    kwargs = {}
    if 'preserve_literal_na_strings' in assignment:
        if assignment['preserve_literal_na_strings'] is not True:
            raise ValueError('Literal metadata token preservation requires explicit true')
        # Source labels such as treatment "None" can be meaningful. Reviewed
        # opt-in preserves them while empty CSV fields remain true missingness.
        kwargs = dict(keep_default_na=False, na_values=[''])
    frame = pd.read_csv(io.BytesIO(raw), index_col=0, dtype=str, **kwargs)
    prefix = assignment['barcode_prefix']
    if not frame.index.is_unique or not frame.index.str.startswith(prefix).all():
        raise ValueError('Duplicate or unexpected source metadata barcodes')
    frame.index = frame.index.str.slice(len(prefix))
    if 'barcode_suffix_replace' in assignment:
        old, new = assignment['barcode_suffix_replace']
        if not old or not frame.index.str.endswith(old).all(): raise ValueError('Unexpected source barcode suffix')
        frame.index = frame.index.str.slice(stop=-len(old)) + new
        if not frame.index.is_unique: raise ValueError('Barcode replacement creates duplicates')
    keys = pd.Index(barcodes.astype(str))
    if not keys.is_unique or not frame.index.isin(keys).all():
        raise ValueError('Metadata barcodes absent from matrix or matrix duplicates')
    donors = frame[assignment['donor_column']]
    if donors.isna().any() or donors.str.strip().eq('').any():
        raise ValueError('Missing donor assignment')
    if set(donors) != set(assignment['expected_donors']):
        raise ValueError('Source donor set changed')
    excluded = donors.isin(assignment.get('excluded_donors', []))
    eligible = frame.loc[~excluded]
    keep = keys.isin(eligible.index)
    audit = dict(matrix_cells=len(keys), source_assigned_cells=len(frame),
                 unassigned_cells=int((~keys.isin(frame.index)).sum()),
                 excluded_donor_cells=int(excluded.sum()), retained_cells=int(keep.sum()))
    return keep, eligible.reindex(keys[keep]), audit


def convert_geo_tar(archive, destination, cfg, metadata_overrides=None):
    """A reviewed explicit member triplet and donor mapping is required for each sample."""
    objects = []; reference_ids = None; assignment_audits = {}; preconversion_audits = {}
    with tarfile.open(archive) as tar:
        for spec in cfg['geo_samples']:
            nested_content = {}
            if 'nested_archive_member' in spec:
                outer = tar.getmember(spec['nested_archive_member'])
                if not outer.isfile() or outer.size > 5_000_000_000:
                    raise ValueError('Invalid nested count archive')
                with tarfile.open(fileobj=io.BytesIO(tar.extractfile(outer).read()), mode='r:*') as inner:
                    for key in ['genes_member', 'barcodes_member', 'matrix_member']:
                        name = spec[key]
                        matches = [m for m in inner.getmembers() if m.name == name]
                        if len(matches) != 1 or not matches[0].isfile() or matches[0].size > 5_000_000_000:
                            raise ValueError('Invalid or ambiguous nested matrix member')
                        nested_content[name] = inner.extractfile(matches[0]).read()
            def content(name):
                if metadata_overrides is not None and name in metadata_overrides:
                    if name != spec.get('cell_assignment', {}).get('member'):
                        raise ValueError('Overrides may supply cell metadata only')
                    raw = metadata_overrides[name]
                    return gzip.decompress(raw) if raw.startswith(b'\x1f\x8b') else raw
                if nested_content:
                    raw = nested_content[name]
                    return gzip.decompress(raw) if raw.startswith(b'\x1f\x8b') else raw
                member = tar.getmember(name)
                if not member.isfile(): raise ValueError('Expected a regular archive member')
                if member.size > 5_000_000_000: raise ValueError('Unexpectedly large matrix member')
                handle = tar.extractfile(member)
                raw = handle.read()
                if raw.startswith(b'Rar!'):
                    # Some GEO submissions label RAR archives .gz. Read a single
                    # member through Windows bsdtar stdout; never extract its path.
                    work = Path(archive).resolve().parent / ('inner_rar_' + uuid.uuid4().hex)
                    work.mkdir()
                    try:
                        rar = Path(work) / 'member.rar'; rar.write_bytes(raw)
                        names = subprocess.run(['tar', '-tf', str(rar)], check=True, capture_output=True).stdout.decode().splitlines()
                        if len(names) != 1: raise ValueError('Expected exactly one file in inner RAR archive')
                        return subprocess.run(['tar', '-xOf', str(rar), names[0]], check=True, capture_output=True).stdout
                    finally:
                        if (work / 'member.rar').exists(): (work / 'member.rar').unlink()
                        work.rmdir()
                # GEO can repack decompressed members while retaining their .gz suffix.
                return gzip.decompress(raw) if raw.startswith(b'\x1f\x8b') else raw
            if 'rna_table_member' in spec:
                if 'stream_rna_table' in spec:
                    if spec['stream_rna_table'] is not True or 'rna_table_layout' in spec or not spec['rna_table_member'].endswith('.gz'):
                        raise ValueError('Streamed RNA table requires explicit true, a gzip member and the standard gene-symbol layout')
                    from src.streamed_gene_table import read_gene_table_stream
                    member=tar.getmember(spec['rna_table_member'])
                    if not member.isfile():raise ValueError('RNA table must be a regular archive member')
                    with tar.extractfile(member) as compressed, gzip.GzipFile(fileobj=compressed) as stream:
                        counts,genes,barcodes=read_gene_table_stream(stream,spec.get('rna_table_separator','\t'))
                elif 'rna_table_layout' in spec:
                    if spec['rna_table_layout'] != 'three_gene_identifiers':
                        raise ValueError('Unreviewed gene-table metadata layout')
                    from src.annotated_gene_table import read_table
                    counts, genes, barcodes = read_table(content(spec['rna_table_member']))
                else:
                    counts, genes, barcodes = read_gene_table(content(spec['rna_table_member']), spec.get('rna_table_separator', '\t'))
            elif 'matrix_h5_member' in spec:
                with h5py.File(io.BytesIO(content(spec['matrix_h5_member'])), 'r') as handle:
                    genome = spec.get('legacy_10x_genome')
                    if genome is not None and list(handle.keys()) != [genome]:
                        raise ValueError('Legacy10xgenomegroups differ fromreviewedsource')
                    matrix = handle[genome if genome is not None else 'matrix']
                    if 'gene_feature_selection' in spec:
                        if cfg['study_id']!='GSE235512':raise ValueError('Gene/TE source selection requires its original study review')
                        from src.gene_te_h5 import read_matrix
                        counts,genes,barcodes,selection_review,source_libraries=read_matrix(matrix)
                        if selection_review!=spec['gene_feature_selection']:raise ValueError('Original gene/TE feature-selection audit changed')
                    else:
                        counts = sparse.csc_matrix((matrix['data'][()], matrix['indices'][()], matrix['indptr'][()]), shape=tuple(matrix['shape'][()])).T.tocsr()
                        def strings(dataset): return [x.decode() if isinstance(x, bytes) else str(x) for x in dataset[()]]
                        if genome is not None:
                            genes = pd.DataFrame({0: strings(matrix['genes']), 1: strings(matrix['gene_names'])})
                            features = {}
                        else:
                            features = matrix['features']
                            genes = pd.DataFrame({0: strings(features['id']), 1: strings(features['name'])})
                        barcodes = pd.Series(strings(matrix['barcodes']))
                        if 'feature_type' in features:
                            rna = np.array(strings(features['feature_type'])) == 'Gene Expression'
                            genes = genes.loc[rna].reset_index(drop=True); counts = counts[:, rna]
            else:
                # Literal identifiers such as influenza neuraminidase "NA"
                # are coordinates, not pandas missing-value markers.
                genes = pd.read_csv(io.BytesIO(content(spec['genes_member'])), sep='\t', header=None, dtype=str, keep_default_na=False)
                if genes.shape[1] == 1:
                    if spec.get('gene_identifier_format') != 'single_column_gene_symbols':
                        raise ValueError('Single-column genes require explicit gene-symbol review')
                    if genes[0].isna().any() or genes[0].str.strip().eq('').any():
                        raise ValueError('Missing source gene symbol')
                    genes[1] = genes[0]
                barcodes = pd.read_csv(io.BytesIO(content(spec['barcodes_member'])), sep='\t', header=None, dtype=str)[0]
                matrix_raw = content(spec['matrix_member'])
                if spec.get('matrix_missing_dimensions') is True:
                    matrix_raw = restore_missing_dimensions(matrix_raw, len(genes), len(barcodes))
                counts = sparse.csr_matrix(scipy_io.mmread(io.BytesIO(matrix_raw), spmatrix=True).T)
                if 'host_gene_selection' in spec:
                    from src.human_ensembl_rna import select
                    counts, genes, source_libraries = select(counts, genes, spec['host_gene_selection'])
                if spec.get('non_rna_feature_ids'):
                    from src.original_triplet_counts import explicit_rna_mask
                    if counts.shape != (len(barcodes),len(genes)) or not np.isfinite(counts.data).all() or (counts.data < 0).any() or not np.equal(counts.data,np.floor(counts.data)).all():
                        raise ValueError('Invalid original counts before explicit modality selection')
                    rna = explicit_rna_mask(genes,spec['non_rna_feature_ids'])
                    genes = genes.loc[rna].reset_index(drop=True); counts = counts[:,rna]
                if genes.shape[1] >= 3:
                    rna = genes[2].eq('Gene Expression').to_numpy()
                    if not rna.any(): raise ValueError('No Gene Expression features in source matrix')
                    genes = genes.loc[rna].reset_index(drop=True); counts = counts[:, rna]
            if counts.shape != (len(barcodes), len(genes)):
                raise ValueError('Matrix dimensions do not match genes and barcodes')
            if genes[0].isna().any() or genes[0].str.strip().eq('').any(): raise ValueError('Missing source gene ID')
            if not genes[0].is_unique: raise ValueError('Duplicate source gene IDs require explicit review')
            ids = genes[0].to_numpy()
            if reference_ids is not None and not np.array_equal(ids, reference_ids):
                raise ValueError('Differing source gene universes require per-sample availability handling')
            reference_ids = ids
            var = pd.DataFrame({'gene_symbols': genes[1].to_numpy(), 'feature_name': genes[1].to_numpy()}, index=ids)
            if spec.get('rna_table_layout') == 'three_gene_identifiers':
                var['source_gene_secondary_id'] = genes['source_gene_secondary_id'].to_numpy()
            if 'gene_symbol_aliases' in spec:
                aliases = spec['gene_symbol_aliases']
                if 'gene_symbol_prefix' in spec or not isinstance(aliases, dict) or not aliases:
                    raise ValueError('Explicit nonempty symbol aliases cannot be combined with prefix mapping')
                if not set(aliases).issubset(set(genes[1])) or any(not isinstance(v,str) or not v.strip() for v in aliases.values()):
                    raise ValueError('Invalid or absent source symbol aliases')
                var['source_gene_symbols'] = genes[1].to_numpy()
                var['gene_symbols'] = genes[1].map(lambda value: aliases.get(value,value)).to_numpy()
                var['feature_name'] = var['gene_symbols'].to_numpy()
            if 'gene_symbol_prefix' in spec:
                prefix = spec['gene_symbol_prefix']
                symbols = genes[1].astype(str)
                if not isinstance(prefix, str) or not prefix or not symbols.str.startswith(prefix).all():
                    raise ValueError('Source gene symbols do not match reviewed prefix')
                symbols = symbols.str.slice(len(prefix))
                if symbols.eq('').any():
                    raise ValueError('Empty gene symbol after reviewed prefix removal')
                var['source_gene_symbols'] = genes[1].to_numpy()
                var['gene_symbols'] = symbols.to_numpy()
                var['feature_name'] = symbols.to_numpy()
            if spec.get('exclude_ambiguous_barcodes',False):
                if 'cell_assignment' in spec or 'preconversion_min_rna_counts' in cfg:
                    raise ValueError('Ambiguous-barcode exclusion needs its own reviewed unpooled source selection')
                if not np.isfinite(counts.data).all() or (counts.data<0).any() or not np.equal(counts.data,np.floor(counts.data)).all():
                    raise ValueError('Invalid original RNA counts before ambiguous-barcode exclusion')
                ambiguous=barcodes.duplicated(keep=False).to_numpy()
                review=spec.get('source_count_review',{})
                if len(barcodes)!=review.get('source_barcode_rows') or int(ambiguous.sum())!=review.get('excluded_ambiguous_barcode_rows'):
                    raise ValueError('Audited source barcode exclusion changed')
                if int(counts.sum())!=review.get('source_total_RNA_counts') or int(counts[ambiguous].sum())!=review.get('excluded_ambiguous_RNA_counts'):
                    raise ValueError('Audited original/excluded RNA count totals changed')
                preconversion_audits[spec['accession']]=dict(source_barcode_rows=len(barcodes),
                    excluded_ambiguous_barcode_rows=int(ambiguous.sum()),kept_barcode_rows=int((~ambiguous).sum()),
                    source_total_RNA_counts=int(counts.sum()),excluded_ambiguous_RNA_counts=int(counts[ambiguous].sum()),
                    policy='Exclude every column with a repeated original barcode; never rename or choose one copy')
                counts=counts[~ambiguous];barcodes=barcodes.iloc[np.flatnonzero(~ambiguous)].reset_index(drop=True)
            if 'preconversion_min_rna_counts' in cfg:
                threshold = cfg['preconversion_min_rna_counts']
                if 'cell_assignment' in spec or threshold != cfg['qc']['min_counts'] or threshold <= 0:
                    raise ValueError('Early RNA filtering requires unpooled records and the unchanged configured QC threshold')
                if not np.isfinite(counts.data).all() or (counts.data < 0).any() or not np.equal(counts.data, np.floor(counts.data)).all():
                    raise ValueError('Invalid original RNA counts before early QC')
                totals = np.asarray(counts.sum(axis=1)).ravel()
                keep = totals >= threshold
                preconversion_audits[spec['accession']] = dict(source_barcode_rows=len(barcodes),
                    RNA_positive_barcode_rows=int((totals > 0).sum()),
                    excluded_below_min_counts=int((~keep).sum()), kept_barcode_rows=int(keep.sum()),
                    min_rna_counts=threshold)
                counts = counts[keep]
                barcodes = barcodes.iloc[np.flatnonzero(keep)].reset_index(drop=True)
            obs = pd.DataFrame(index=spec['accession'] + ':' + barcodes.astype(str))
            obs.index.name = 'cell_id'
            obs['source_sample_id'] = spec['accession']; obs['donor_id'] = spec['donor_id']
            if 'gene_feature_selection' in spec or 'host_gene_selection' in spec:
                selected_libraries=source_libraries.loc[keep].reset_index(drop=True) if 'preconversion_min_rna_counts' in cfg else source_libraries
                for name in selected_libraries:obs[name]=selected_libraries[name].to_numpy()
            obs['tissue'] = 'blood'; obs['source_specimen'] = 'PBMC'; obs['is_primary_data'] = True
            for key, value in spec['metadata'].items():
                # A JSON null is true missingness. Numeric values from other
                # captures must not turn this into an unserializable object.
                obs[key] = np.nan if value is None else value
            if 'cell_assignment' in spec:
                assignment = spec['cell_assignment']
                keep, assigned, audit = align_pool_metadata(content(assignment['member']), barcodes, assignment)
                counts = counts[keep]; obs = obs.iloc[np.flatnonzero(keep)].copy()
                obs['donor_id'] = assigned[assignment['donor_column']].to_numpy()
                obs['source_sample_id'] = obs['donor_id'].to_numpy()
                obs['source_library_id'] = spec['accession']
                for key in assigned.columns:
                    # CSV assignment fields are strings or missing. Categorical
                    # codes preserve true missingness even for all-empty columns
                    # and avoid HDF5 trying to serialize float NaNs as strings.
                    obs['source_cell_' + key.replace('.', '_')] = pd.Categorical(assigned[key])
                assignment_audits[spec['accession']] = audit
            obj = ad.AnnData(counts, obs=obs, var=var)
            objects.append(obj)
    combined = ad.concat(objects, join='inner', merge='same')
    # Keep the verified inner gene axis while retaining fields supplied for
    # only some source captures. Missing metadata is not a reason to drop a
    # whole clinical/provenance column from every capture.
    complete_obs = pd.concat([obj.obs for obj in objects],axis=0,join='outer')
    for column in complete_obs.columns:
        if complete_obs[column].dtype == object:
            kind = pd.api.types.infer_dtype(complete_obs[column].dropna(),skipna=True)
            if kind in {'string','unicode','boolean'}:
                complete_obs[column] = pd.Categorical(complete_obs[column])
    combined.obs = complete_obs.loc[combined.obs_names].copy()
    # These cohorts' X is already the original count matrix. Keeping a second
    # identical copy in raw.X would add large temporary disk use without data.
    matrix_bytes = combined.X.data.nbytes + combined.X.indices.nbytes + combined.X.indptr.nbytes
    # Large GEO assemblies already contain original counts in X. A second raw.X
    # payload can exceed the project budget without adding any information.
    duplicate_raw = (cfg.get('dataset_id') not in {'GSE309595_genes001', 'GSE318371'}
                     and matrix_bytes < 1_000_000_000)
    if duplicate_raw:
        combined.raw = combined.copy()
    else:
        combined.uns['count_carrier_X_original'] = True
    combined.uns['input_note'] = ('Source raw count matrices in X; no pre-normalization. '
        'All source metadata preserved in configuration. ' +
        ('raw.X duplicates X.' if duplicate_raw else 'raw.X omitted because it would duplicate X.'))
    if assignment_audits:
        combined.uns['source_assignment_audit'] = assignment_audits
    if preconversion_audits:
        combined.uns['source_preconversion_audit_json'] = json.dumps(preconversion_audits)
    from scripts.download_pbmc_asset import footprint, LIMIT
    # Reserve the uncompressed sparse arrays plus metadata and HDF5 overhead
    # before writing the converted temporary input.
    # With one count matrix, H5AD writes one compressed sparse payload. This
    # large in-memory sparse arrays derive from compressed source counts.
    # Use a bounded prewrite estimate and verify the actual
    # project footprint before publishing the converted carrier.
    target = Path(destination)
    if not duplicate_raw and target.exists():
        # A prior interrupted replay may have left its converted count carrier.
        # Rebuild from the freshly checksum-verified GEO assets above and compare
        # every cell and gene before reusing it. Never overwrite an unknown file.
        existing = ad.read_h5ad(target, backed='r')
        try:
            if (existing.shape != combined.shape or
                    (existing.raw is None and not bool(existing.uns.get('count_carrier_X_original', False))) or
                    not np.array_equal(existing.obs_names.astype(str), combined.obs_names.astype(str)) or
                    not np.array_equal(existing.var_names.astype(str), combined.var_names.astype(str))):
                raise ValueError('Existing converted carrier axes or count location differ')
            if (existing.raw is not None and
                    not np.array_equal(existing.raw.var_names.astype(str), combined.var_names.astype(str))):
                raise ValueError('Existing converted carrier raw gene axis differs')
            for field in ('source_sample_id', 'donor_id'):
                if (field not in existing.obs or field not in combined.obs or
                        not np.array_equal(existing.obs[field].astype(str),
                                           combined.obs[field].astype(str))):
                    raise ValueError(f'Existing converted carrier {field} differs')
            for start in range(0, combined.n_obs, 256):
                stop = min(start + 256, combined.n_obs)
                current = sparse.csr_matrix(existing.X[start:stop])
                rebuilt = sparse.csr_matrix(combined.X[start:stop])
                if (current - rebuilt).nnz:
                    raise ValueError(f'Existing converted carrier counts differ at row {start}')
                if existing.raw is not None:
                    raw = sparse.csr_matrix(existing.raw.X[start:stop])
                    if (raw - rebuilt).nnz:
                        raise ValueError(f'Existing converted carrier raw counts differ at row {start}')
        finally:
            existing.file.close()
        return target
    obs_bytes = int(combined.obs.memory_usage(deep=True).sum())
    estimated = (3 if duplicate_raw else 0.6) * matrix_bytes + 2 * obs_bytes + (100_000_000 if duplicate_raw else 2_000_000_000)
    used = footprint()
    if used + estimated >= LIMIT:
        raise RuntimeError(f'Converted input would exceed the project disk budget: used={used}, '
                           f'estimated={int(estimated)}, matrix={matrix_bytes}, obs={obs_bytes}, limit={LIMIT}')
    if duplicate_raw:
        combined.write_h5ad(destination, compression='gzip')
    else:
        temporary = target.with_suffix('.building.h5ad')
        if temporary.exists() or target.exists():
            raise RuntimeError('Inspect existing converted carrier before retrying')
        try:
            combined.write_h5ad(temporary, compression='gzip')
            actual = footprint()
            if actual >= LIMIT:
                raise RuntimeError(f'Converted input exceeded project budget after write: {actual} >= {LIMIT}')
            temporary.replace(target)
        except Exception:
            if temporary.exists():
                temporary.unlink()
            raise
    return Path(destination)
