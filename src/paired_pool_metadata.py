"""Reconcile the two published GSE194315 cell assignment/QC tables."""
import copy
import json
from pathlib import Path
import pandas as pd
from src.geo_inputs import convert_geo_tar
from scripts.download_pbmc_asset import ROOT, download


def select_cells(first, second):
    a = first.set_index('CellName'); b = second.set_index('CellName')
    if not a.index.is_unique or not b.index.is_unique or set(a.index) != set(b.index):
        raise ValueError('Source cell identifiers differ or duplicate')
    b = b.reindex(a.index)
    for col in ['Sample', 'Run', 'Status', 'DemuxletDropletType']:
        if not a[col].equals(b[col]): raise ValueError('Conflicting source field: ' + col)
    disagreement = b.Subject.notna() & a.Subject.ne(b.Subject)
    if disagreement.any(): raise ValueError('Conflicting nonmissing source donor assignments')
    for frame in [a, b]:
        if not frame.IncludedInStudy.isin(['TRUE', 'FALSE']).all():
            raise ValueError('Unrecognized source QC flags')
    keep = a.DemuxletDropletType.eq('SNG') & (a.IncludedInStudy.eq('TRUE') | b.IncludedInStudy.eq('TRUE')) & a.Subject.ne('X01')
    result = a.loc[keep].copy()
    if result.Subject.isna().any(): raise ValueError('Missing selected donor')
    result['IncludedInASStudy'] = a.loc[keep, 'IncludedInStudy']
    result['IncludedInPSAStudy'] = b.loc[keep, 'IncludedInStudy']
    return result


def prepare(cfg, destination):
    paths = {}; manifests = {}
    for role, asset in cfg['source_assets'].items():
        name = asset['filename']
        if Path(name).name != name: raise ValueError('Unsafe filename')
        paths[role] = download(asset['url'], ROOT / 'data/downloads' / name, asset['expected_bytes'])
        manifests[role] = json.loads(paths[role].with_suffix(paths[role].suffix + '.manifest.json').read_text())
    selected = select_cells(pd.read_csv(paths['as_metadata'], sep='\t', dtype=str),
                            pd.read_csv(paths['psa_metadata'], sep='\t', dtype=str))
    if set(selected.Subject) != set(cfg['expected_donors']): raise ValueError('Donor set changed')
    if set(selected.Sample) != {s['accession'] for s in cfg['geo_samples']}: raise ValueError('Library set changed')
    derived = copy.deepcopy(cfg); overrides = {}
    for spec in derived['geo_samples']:
        library = spec['accession']; frame = selected.loc[selected.Sample.eq(library)].copy()
        # Source metadata omits the standard 10x -1 suffix; matrix barcodes include it.
        if not frame.index.str.startswith(library + '_').all(): raise ValueError('Unexpected source barcode prefix')
        bare = frame.index.str.slice(len(library) + 1)
        if not bare.str.fullmatch('[ACGT]{16}').all(): raise ValueError('Unexpected source barcode format')
        frame.index = library + '_' + bare + '-1'
        name = library + '_reviewed_assignments.csv'
        overrides[name] = frame.to_csv().encode()
        spec['cell_assignment'] = dict(member=name, barcode_prefix=library + '_', donor_column='Subject',
            expected_donors=sorted(frame.Subject.unique()))
    convert_geo_tar(paths['counts'], destination, derived, metadata_overrides=overrides)
    return dict(kind='GSE194315_paired_metadata', assets=manifests, selected_cells=len(selected),
        selected_donors=selected.Subject.nunique(), source_qc='Singlets passing either author QC table; X01 excluded',
        barcode_mapping='Exact library prefix and 16nt barcode, append literal -1 to match deposited matrix'), [paths['counts']]
