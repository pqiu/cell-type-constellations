"""Read an explicitly reviewed Seurat RNA count slot and source cell metadata."""
import json
from pathlib import Path
import anndata as ad
import pandas as pd
import rdata
from src.rdata_counts import sparse_counts
from src.source_selection import source_exclusion_mask
from scripts.download_pbmc_asset import ROOT, LIMIT, footprint, download


def read_count_slots(path, assay, compression_layers=None):
    """Convert only named counts and metadata; avoid unrelated R command closures."""
    if compression_layers is None:
        parsed = rdata.parser.parse_file(path)
    else:
        from src.rds_streaming import parse_file
        parsed = parse_file(path, compression_layers)
    def attributes(obj):
        result = {}; node = obj.attributes; seen = set()
        while node is not None and node.info.type != rdata.parser.RObjectType.NILVALUE:
            if id(node) in seen or node.info.type != rdata.parser.RObjectType.LIST:
                raise ValueError('Invalid R attribute pairlist')
            seen.add(id(node))
            name = rdata.conversion.convert(node.tag)
            if name in result: raise ValueError('Duplicate R attribute')
            result[name] = node.value[0]; node = node.value[1]
        return result
    root = attributes(parsed.object)
    if 'Seurat' not in rdata.conversion.convert(root['class']):
        raise ValueError('Expected Seurat object')
    assays = root['assays']
    names = list(rdata.conversion.convert(attributes(assays)['names']))
    if len(names) != len(set(names)) or len(names) != len(assays.value) or assay not in names:
        raise ValueError('Invalid assay names')
    selected = attributes(assays.value[names.index(assay)])
    counts = rdata.conversion.convert(selected['counts'])
    metadata = rdata.conversion.convert(root['meta.data'])
    return counts, metadata


def prepare(cfg, destination):
    asset = cfg['source_asset']; name = asset['filename']
    if Path(name).name != name: raise ValueError('Unsafe source name')
    path = download(asset['url'], ROOT / 'data/downloads' / name, asset['expected_bytes'])
    manifest = json.loads(path.with_suffix(path.suffix+'.manifest.json').read_text())
    if cfg.get('selective_rds_slots', False):
        raw_counts, metadata = read_count_slots(path, cfg['source_assay'], cfg.get('rds_compression_layers'))
    else:
        source = rdata.read_rds(path)
        raw_counts = source.assays[cfg['source_assay']].counts
        metadata = getattr(source, 'meta.data')
    counts, genes, cells = sparse_counts(raw_counts)
    obs = metadata.copy()
    if not obs.index.is_unique or set(obs.index) != set(cells): raise ValueError('Source cell identities differ')
    obs = obs.loc[cells]
    excluded = source_exclusion_mask(obs, cfg.get('source_cell_exclusions', {}))
    obs = obs.loc[~excluded].copy(); counts = counts[~excluded].copy()
    for column in obs:
        if obs[column].dtype == object or isinstance(obs[column].dtype, pd.StringDtype):
            obs[column] = obs[column].fillna('not_reported').astype(str)
    if footprint() + 2*(counts.data.nbytes+counts.indices.nbytes+counts.indptr.nbytes) + 2_000_000_000 >= LIMIT:
        raise RuntimeError('Seurat count conversion exceeds disk budget')
    var = pd.DataFrame({'gene_symbols': genes.to_numpy()}, index=genes)
    ad.AnnData(X=counts, obs=obs, var=var).write_h5ad(destination, compression='gzip')
    return dict(kind='reviewed_seurat_RNA_counts', source=manifest, reader_version=rdata.__version__,
        source_cells=len(cells), excluded_source_cells=int(excluded.sum()), selected_cells=len(obs),
        conversion='RNA counts slot only; source metadata barcode aligned; no normalization or embedding used.'), [path]
