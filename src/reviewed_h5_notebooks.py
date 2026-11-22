"""Validate original modern10x RNA H5 assets for explicitly reviewed specimens."""
import json
import h5py
import pandas as pd
from scripts.download_pbmc_asset import ROOT, download
from src.lyme_vireo import audit_h5
from src.reviewed_triplet_notebooks import export_groups


def prepare(accession, specimens, title, evidence, publications):
    if not specimens or len({s['accession'] for s in specimens}) != len(specimens):
        raise ValueError('Unique reviewed original records required')
    groups = {}; audits = {}; seen = set()
    for spec in specimens:
        row = spec['source_metadata']; gsm = spec['accession']
        if row['geo_accession'] != [gsm] or row['organism_ch1'] != ['Homo sapiens']:
            raise ValueError('Original source organism/record mismatch')
        urls = [u.replace('ftp://', 'https://') for k, vs in row.items() if k.startswith('supplementary_file') for u in vs]
        hits = [u for u in urls if u.endswith('filtered_feature_bc_matrix.h5')]
        if len(hits) != 1: raise ValueError('Ambiguous original RNA H5')
        url = hits[0]; name = url.rsplit('/', 1)[1]
        if not name.startswith(gsm + '_'): raise ValueError('Original source asset mismatch')
        path = download(url, ROOT / f'data/downloads/{accession}_assets' / name)
        with h5py.File(path, 'r') as f:
            barcodes = f['matrix/barcodes'].asstr()[()]
            types = f['matrix/features/feature_type'].asstr()[()]
        review = audit_h5(path, pd.DataFrame(index=barcodes), seen)
        review['source_feature_types'] = pd.Series(types).value_counts().to_dict()
        spec.update(matrix_h5_member=name, assets=[dict(url=url, filename=name)], source_count_review=review)
        audits[gsm] = review; groups.setdefault(review['gene_coordinate_sha256'], []).append(spec)
        print(dict(specimen=gsm, donor=spec['donor_id'], source_barcodes=len(barcodes)), flush=True)
    export_groups(accession, specimens, title, evidence, publications, groups, audits)

