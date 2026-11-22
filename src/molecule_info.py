"""Recover RNA UMI counts from reviewed per-sample Cell Ranger molecule files."""
import json
import h5py
import numpy as np
import pandas as pd
from scipy import sparse


def read_counts(path, chunk_size=2_000_000, reviewed_source_mode=None):
    with h5py.File(path) as f:
        metrics = json.loads(f['metrics_json'][()].decode())
        exonic_count = reviewed_source_mode == 'cellranger6_count_exonic'
        if reviewed_source_mode not in (None, 'cellranger6_count_exonic'):
            raise ValueError('Unknown reviewed molecule source mode')
        if exonic_count:
            if (metrics.get('molecule_info_type') != 'count'
                    or metrics.get('cellranger_version') != 'cellranger-6.0.2'
                    or metrics.get('analysis_parameters', {}).get('include_introns') is not False
                    or int(f.attrs.get('file_version', -1)) != 5
                    or f.attrs.get('filetype') != 'molecule'
                    or metrics.get('is_aggregated', False) is not False
                    or 'umi_type' not in f):
                raise ValueError('Expected explicitly reviewed Cell Ranger6.0.2 exonic count output')
        else:
            if metrics.get('molecule_info_type') != 'per_sample':
                raise ValueError('Reviewed per-sample molecule output required')
            if metrics['analysis_parameters'].get('include_introns') is not True:
                raise ValueError('Only reviewed intron-inclusive inputs supported')
        libraries = json.loads(f['library_info'].asstr()[0])
        rna = [i for i, x in enumerate(libraries) if x['library_type'] == 'Gene Expression']
        if len(rna) != 1:
            raise ValueError('Exactly one RNA library required')
        features = f['features/feature_type'].asstr()[:] == 'Gene Expression'
        ids = f['features/id'].asstr()[:][features]
        symbols = f['features/name'].asstr()[:][features]
        if len(set(ids)) != len(ids):
            raise ValueError('Duplicate RNA gene coordinates')
        if set(f['features/genome'].asstr()[:][features]) != {'GRCh38'}:
            raise ValueError('Expected reviewed GRCh38 RNA reference')
        pf = f['barcode_info/pass_filter'][:]
        pf = pf[pf[:, 1] == rna[0]]
        if len(pf) == 0 or len(set(pf[:, 0])) != len(pf) or set(pf[:, 2]) != {0}:
            raise ValueError('Ambiguous filtered RNA cell identities')
        allowed = np.sort(pf[:, 0])
        barcodes = f['barcodes'].asstr()[:][allowed]
        if len(set(barcodes)) != len(barcodes):
            raise ValueError('Duplicate source barcodes')
        barcode_ok = np.zeros(len(f['barcodes']), dtype=bool)
        barcode_ok[allowed] = True
        parts = []; supporting_reads = 0; intronic_filtered_RNA_molecules = 0
        for start in range(0, len(f['umi']), chunk_size):
            sl = slice(start, start+chunk_size)
            b = f['barcode_idx'][sl]; g = f['feature_idx'][sl]
            counts = f['count'][sl]; lib = f['library_idx'][sl]
            keep = (lib == rna[0]) & barcode_ok[b] & features[g] & (counts > 0)
            if exonic_count:
                types = f['umi_type'][sl]
                if not np.isin(types[lib == rna[0]], [0, 1]).all():
                    raise ValueError('Unknown RNA molecule type bits')
                intronic_filtered_RNA_molecules += int(np.count_nonzero(keep & (types == 0)))
                keep &= types == 1
            if not np.all(f['gem_group'][sl][keep] == 1):
                raise ValueError('Multiple GEM groups require explicit barcode mapping')
            keys = np.empty(int(keep.sum()), dtype=[('barcode','<u4'), ('gene','<u4'), ('umi','<u8')])
            keys['barcode'] = b[keep]; keys['gene'] = g[keep]; keys['umi'] = f['umi'][sl][keep]
            parts.append(keys); supporting_reads += int(counts[keep].sum())
        keys = np.concatenate(parts)
        unique = np.unique(keys)
        if len(keys) != len(unique):
            raise ValueError('Duplicate RNA molecule identities require review')
        gene_map = np.full(len(features), -1, dtype=np.int32)
        gene_map[features] = np.arange(len(ids))
        x = sparse.coo_matrix((np.ones(len(keys), dtype=np.int64),
            (np.searchsorted(allowed, keys['barcode']), gene_map[keys['gene']])),
            shape=(len(allowed), len(ids))).tocsr()
        if int(x.sum()) != len(keys):
            raise ValueError('RNA UMI totals do not reconcile')
        var = pd.DataFrame({'gene_symbols': symbols, 'feature_type': 'Gene Expression'}, index=ids)
        report = dict(filtered_cells=len(allowed), rna_features=len(ids), rna_umi_counts=len(keys),
                      supporting_reads=supporting_reads, source_metrics=metrics,
                      count_rule=('One positive-read distinct barcode-gene-UMI record; RNA library and filtered cells only; '
                                  + ('source exonic umi_type1 only, introns excluded.' if exonic_count else 'source introns included.')))
        if exonic_count:
            report['reviewed_source_mode'] = reviewed_source_mode
            report['intronic_positive_RNA_molecules_in_filtered_cells_excluded'] = intronic_filtered_RNA_molecules
        return x, var, pd.Index([b+'-1' for b in barcodes]), report


def prepare(cfg, destination):
    import anndata as ad
    from pathlib import Path
    from scripts.download_pbmc_asset import ROOT, download, footprint, LIMIT
    asset = cfg['source_asset']
    name = asset['filename']
    if Path(name).name != name:
        raise ValueError('Source asset must be a basename')
    subdirectory = cfg.get('source_asset_subdirectory', cfg['study_id'] + '_assets')
    if subdirectory and (Path(subdirectory).name != subdirectory or subdirectory in ['.', '..']):
        raise ValueError('Source asset directory must stay inside downloads')
    path = download(asset['url'], ROOT/'data/downloads'/subdirectory/name, asset.get('expected_bytes'))
    manifest = json.loads(path.with_suffix(path.suffix+'.manifest.json').read_text())
    x, var, cells, report = read_counts(path, reviewed_source_mode=cfg.get('molecule_source_mode'))
    if footprint() + x.data.nbytes + x.indices.nbytes + x.indptr.nbytes + 2_000_000_000 >= LIMIT:
        raise RuntimeError('Converted counts would exceed disk limit')
    ad.AnnData(X=x, obs=pd.DataFrame(index=cells), var=var).write_h5ad(destination, compression='gzip')
    return dict(kind='reviewed_per_sample_RNA_molecules', asset=manifest, conversion=report), [path]
