"""Conservative sample-tag calls; probabilities are not calibrated accuracy."""
import anndata as ad
import numpy as np
import pandas as pd
from scanpy.external.pp import hashsolo

PRIORS = [(0.01, 0.8, 0.19), (0.05, 0.75, 0.20)]


def select(counts, mapping, minimum_probability=0.95):
    if not counts.index.is_unique or not all(isinstance(b, str) and b.strip() for b in counts.index) or len(counts.columns) < 2 or set(counts.columns) != set(mapping):
        raise ValueError('Unique barcodes and exact source tag mapping required')
    values = counts.to_numpy(dtype=float)
    if not np.isfinite(values).all() or (values < 0).any() or not np.equal(values, np.floor(values)).all():
        raise ValueError('Sample tags must be native nonnegative integer counts')
    if not 0.5 <= minimum_probability <= 1: raise ValueError('Invalid posterior threshold')
    # HashSolo's seeded fits and tie handling depend on column order. Canonical
    # tag-ID order makes the same native counts reproduce the same identities.
    counts = counts.reindex(columns=sorted(mapping))
    outputs = []
    # Use HashSolo's n_tags-2 background setting where possible. Two-tag pools
    # explicitly retain one low tag because the automatic setting would be zero.
    noise_barcodes = max(1, len(counts.columns) - 2)
    for priors in PRIORS:
        np.random.seed(0)
        obj = ad.AnnData(obs=counts.copy())
        hashsolo(obj, list(counts.columns), priors=priors, number_of_noise_barcodes=noise_barcodes)
        outputs.append(obj.obs)
    first, second = outputs
    probability_columns = ['negative_hypothesis_probability', 'singlet_hypothesis_probability', 'doublet_hypothesis_probability']
    for result in outputs:
        probabilities = result[probability_columns].to_numpy()
        if not np.isfinite(probabilities).all() or (probabilities < 0).any() or (probabilities > 1 + 1e-12).any():
            raise ValueError('Invalid assignment posterior')
        if not np.allclose(probabilities.sum(axis=1), 1, atol=1e-10):
            raise ValueError('Assignment probabilities do not sum to one')
        if not set(result.Classification).issubset(set(mapping) | {'Negative', 'Doublet'}):
            raise ValueError('Unexpected tag call')
    retained = first.Classification.eq(second.Classification) & first.Classification.isin(mapping)
    retained &= first.singlet_hypothesis_probability.ge(minimum_probability)
    retained &= second.singlet_hypothesis_probability.ge(minimum_probability)
    calls = first.copy()
    calls['sensitivity_classification'] = second.Classification
    calls['sensitivity_singlet_probability'] = second.singlet_hypothesis_probability
    calls['assignment_retained'] = retained
    calls['specimen_id'] = calls.Classification.map(mapping).where(retained)
    calls['assignment_method'] = 'HashSolo; two-prior agreement; model posterior not calibrated correctness'
    return calls, calls.loc[retained].copy()
