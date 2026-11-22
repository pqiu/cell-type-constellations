"""Read an explicitly paired count matrix, gene table and cell metadata table.

Matrix row/column order must be documented by the original deposition. Numeric
aggregation ordinals never supply donor identities. Derived axis files merely
adapt the existing strict integer-count converter and are temporary inputs.
"""
import json
from pathlib import Path

import pandas as pd

from scripts.download_pbmc_asset import ROOT, download
from src.geo_mtx_table import convert


def prepare(cfg, destination):
    folder = ROOT / 'data/downloads' / cfg['source_asset_subdirectory']
    if Path(cfg['source_asset_subdirectory']).name != cfg['source_asset_subdirectory']:
        raise ValueError('Source directory must be a basename')
    paths = {}; manifests = {}
    if set(cfg['matrix_assets']) != {'counts', 'cells', 'genes'}:
        raise ValueError('Three explicitly paired original assets required')
    if not cfg.get('source_pairing_evidence'):
        raise ValueError('Original row/column pairing evidence is required')
    for role, asset in cfg['matrix_assets'].items():
        name = asset['filename']
        if Path(name).name != name:
            raise ValueError('Original source filename must be a basename')
        p = download(asset['url'], folder / name, asset['expected_bytes'])
        manifest = json.loads(p.with_suffix(p.suffix + '.manifest.json').read_text())
        if manifest['sha256'] != asset['sha256']:
            raise ValueError('Pinned original paired-table source changed')
        paths[role] = p; manifests[role] = manifest
    genes = pd.read_csv(paths['genes'], sep='\t', dtype=str, keep_default_na=False)
    if list(genes) != [cfg['source_gene_column']]:
        raise ValueError('Original gene-table header changed')
    symbols = genes[cfg['source_gene_column']]
    if not symbols.is_unique or symbols.str.strip().eq('').any():
        raise ValueError('Original nonempty unique gene symbols required')
    obs = pd.read_csv(paths['cells'], sep='\t', index_col=0, keep_default_na=False)
    obs.index = obs.index.astype(str)
    if not obs.index.is_unique or obs.index.tolist() != obs[cfg['source_barcode_column']].astype(str).tolist():
        raise ValueError('Original cell-row/barcode pairing changed')
    library_column = cfg['source_library_column']; ordinal_column = cfg['source_ordinal_column']
    if library_column == ordinal_column:
        raise ValueError('Aggregation ordinals cannot supply clinical specimen IDs')
    if set(obs[library_column].astype(str)) != set(cfg['source_all_library_ids']):
        raise ValueError('Original subject-library ID set changed')
    excluded = set(cfg['source_excluded_library_ids'])
    if not excluded.issubset(set(cfg['source_all_library_ids'])):
        raise ValueError('Unknown held source library')
    if 'source_author_aggregation_ordinal' in obs:
        raise ValueError('Refuse original source-field overwrite')
    obs = obs.rename(columns={ordinal_column: 'source_author_aggregation_ordinal'})
    obs['source_subject_library_id'] = obs[library_column].astype(str)
    obs['source_scope'] = obs[library_column].map(lambda key: 'held_identity' if key in excluded else 'PBMC_containing')
    features_path = destination.with_suffix('.derived_features.tsv')
    barcodes_path = destination.with_suffix('.derived_barcodes.tsv')
    metadata_path = destination.with_suffix('.derived_metadata.csv')
    pd.DataFrame({0: symbols, 1: symbols, 2: 'Gene Expression'}).to_csv(features_path, sep='\t', index=False, header=False)
    pd.Series(obs.index).to_csv(barcodes_path, sep='\t', index=False, header=False)
    obs.to_csv(metadata_path)
    conversion = convert(dict(matrix=paths['counts'], features=features_path,
        barcodes=barcodes_path, metadata=metadata_path), destination, cfg)
    return (dict(kind='explicitly_paired_original_matrix_and_tables', assets=manifests,
        source_pairing_evidence=cfg['source_pairing_evidence'], conversion=conversion,
        excluded_identity_libraries=sorted(excluded),
        gene_coordinate_policy='Original symbols retained as IDs and names; no invented Ensembl identifiers.'),
        list(paths.values()) + [features_path, barcodes_path, metadata_path])
