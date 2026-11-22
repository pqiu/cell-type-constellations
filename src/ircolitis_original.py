"""Reassemble the authors' native PBMC RNA partitions using exact source cell keys."""
import json

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from scripts.audit_PBMC_206298_original_counts import expand_source, selected_source
from scripts.download_pbmc_asset import ROOT, LIMIT, download, footprint


def canonical_partitions(parts):
    """Count only exact repeated cell IDs once, preserving every author cluster call."""
    matrices = []; observations = []; axis = None; locations = {}; calls = {}; copies = []
    for label, x, obs, genes in parts:
        if axis is None: axis = genes.copy()
        else: pd.testing.assert_frame_equal(genes, axis)
        if not obs.index.is_unique or x.shape != (len(obs), len(genes)):
            raise ValueError('Unique source cell/gene axes required')
        x = sparse.csr_matrix(x); x.sum_duplicates(); x.sort_indices(); x.eliminate_zeros()
        kept = np.ones(len(obs), dtype=bool)
        retained_position = 0
        for i, cell in enumerate(obs.index):
            cluster = str(obs.iloc[i]['cluster'])
            if cell in locations:
                number, position = locations[cell]
                if (matrices[number][position] != x[i]).nnz:
                    raise ValueError('Shared source cell has different native RNA counts')
                before = observations[number].iloc[position]
                columns = [c for c in obs if c != 'cluster']
                if set(columns) != set(observations[number].columns) - {'cluster'}:
                    raise ValueError('Shared source cell clinical schema changed')
                if any(str(before[c]) != str(obs.iloc[i][c]) for c in columns):
                    raise ValueError('Shared source cell has inconsistent clinical identity')
                kept[i] = False; copies.append(dict(cell=cell, duplicate_partition=label))
            else:
                locations[cell] = (len(matrices), retained_position)
                retained_position += 1
            calls.setdefault(cell, {})[label] = cluster
        matrices.append(x[kept].copy()); observations.append(obs.loc[kept].copy())
    result = pd.concat(observations)
    assert result.index.is_unique
    result = result.drop(columns='cluster')
    result['source_author_cluster_assignments'] = [json.dumps(calls[cell], sort_keys=True) for cell in result.index]
    result['source_author_partitions'] = [';'.join(sorted(calls[cell])) for cell in result.index]
    return sparse.vstack(matrices, format='csr'), result, axis, copies


def prepare(cfg, destination):
    from src.pbmc_pipeline import digest
    def load(pin):
        p = ROOT / pin['path']; assert digest(p) == pin['sha256']
        return json.loads(p.read_text(encoding='utf-8'))
    audit = load(cfg['native_count_audit']); review = load(cfg['source_identity_review'])
    assert review['original_count_audit'] == cfg['native_count_audit']
    assert review['complete'] and audit['all_shared_cell_counts_and_selected_clinical_fields_agree']
    parts = []; paths = []; manifests = []
    for source in audit['source_partitions']:
        asset = source['original_gzip']; path = download(asset['url'], ROOT / asset['path'], asset['bytes'])
        assert digest(path) == asset['sha256']
        expanded = expand_source(path); assert digest(expanded) == source['original_expanded']['sha256']
        x, obs, genes, stored = selected_source(expanded)
        assert len(obs) == source['original_cells'] and stored == source['stored_values_audited']
        parts.append((source['partition'], x, obs, genes)); paths.extend([path, expanded]); manifests.append(asset)
    x, obs, genes, copies = canonical_partitions(parts)
    assert len(obs) == review['unique_original_cells'] == 118817 and len(copies) == 1
    assert obs['sample'].nunique() == obs.donor.nunique() == 21 and obs.followup.eq('Initial').all()
    mapping = load(review['source_initial_sample_crosswalk'])
    mapping = {row['source_patient_id']: row for row in mapping}
    for donor in obs.donor.unique():
        row = mapping[donor]; assert obs.loc[obs.donor.eq(donor), 'sample'].eq(row['source_sample_code']).all()
    obs['source_GEO_accession'] = obs.donor.map({key: value['source_GSM'] for key, value in mapping.items()})
    obs['source_GEO_tissue_as_reported'] = obs.donor.map({key: value['original_GEO_fields'].get('tissue', '') for key, value in mapping.items()})
    obs['source_GEO_disease_state_as_reported'] = obs.donor.map({key: value['original_GEO_fields'].get('disease state', '') for key, value in mapping.items()})
    obs['source_GEO_steroids_at_biopsy_as_reported'] = obs.donor.map({key: value['original_GEO_fields'].get('steroids at_the_time_of_biopsy', '') for key, value in mapping.items()})
    obs['donor_id'] = obs.donor
    obs['tissue'] = 'blood'; obs['is_primary_data'] = True
    obs['source_scope'] = 'Public initial PBMC/blood RNA; retain source sample/capture views separately. Original GEO tissue fields conflict and are preserved.'
    obs['source_day_policy'] = 'Opaque source day/sample clock; not inferred days since therapy or an identifiable collection date.'
    obs['source_capture_policy'] = 'Native batch/capture labels retained separately; no pooling across differently labelled captures or whole-PBMC composition claim.'
    obs['source_split_assignment'] = 'Unassigned. Keep all captures from one donor together; global cross-study aliases unresolved.'
    for column in obs:
        if obs[column].dtype == object: obs[column] = pd.Categorical(obs[column])
    needed = 2 * (x.data.nbytes + x.indices.nbytes + x.indptr.nbytes)
    if footprint() + needed + 2_000_000_000 >= LIMIT: raise RuntimeError('Native RNA assembly exceeds project disk budget')
    ad.AnnData(X=x, obs=obs, var=genes).write_h5ad(destination, compression='gzip')
    return dict(kind='original_ircolitis_RNA_partition_union', source_assets=manifests,
        native_count_audit=cfg['native_count_audit'], source_identity_review=cfg['source_identity_review'],
        exact_source_cell_copies=copies, unique_source_cells=len(obs), native_gene_rows=len(genes),
        source_capture_views=88, count_slot='X', count_transformation='None; exact source rows concatenated, one proved duplicate cell ID counted once.'), paths
