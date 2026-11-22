"""Align source donor assignments to expression barcodes without guessing IDs."""
import numpy as np
import pandas as pd


def align_demuxlet(barcodes, assignments, expected_donors, min_singlet_probability=0.95):
    """Return all expression barcodes and an explicit eligibility/audit flag.

    SNG.POSTERIOR is the probability of being a singlet, not a probability
    that the selected donor identity is correct. Preserve both distinctions.
    Missing assignments remain in the audit and are ineligible.
    """
    barcodes=pd.Index(barcodes)
    if barcodes.isna().any():
        raise ValueError('Expression barcodes must be unique and nonempty')
    barcodes=barcodes.astype(str)
    if barcodes.has_duplicates or barcodes.isna().any() or (barcodes=='').any():
        raise ValueError('Expression barcodes must be unique and nonempty')
    required=['BARCODE','DROPLET.TYPE','SNG.POSTERIOR','SNG.BEST.GUESS']
    if not set(required).issubset(assignments.columns):
        raise ValueError('Missing demuxlet fields')
    if assignments.BARCODE.isna().any() or assignments.BARCODE.astype(str).duplicated().any() or assignments.BARCODE.eq('').any():
        raise ValueError('Assignment barcodes must be unique and nonmissing')
    if not 0 <= min_singlet_probability <= 1:
        raise ValueError('Invalid singlet threshold')
    expected=set(expected_donors)
    if not expected or any(not isinstance(x,str) or not x for x in expected):
        raise ValueError('Expected donor IDs must be explicit nonempty strings')
    table=assignments.copy();table['BARCODE']=table.BARCODE.astype(str)
    prob=pd.to_numeric(table['SNG.POSTERIOR'],errors='raise').to_numpy(dtype=float)
    if not np.isfinite(prob).all() or ((prob<0)|(prob>1)).any():
        raise ValueError('Invalid singlet probability')
    if not table['DROPLET.TYPE'].isin(['SNG','DBL','AMB']).all():
        raise ValueError('Unrecognized source droplet type')
    table['SNG.POSTERIOR']=prob
    aligned=table.set_index('BARCODE').reindex(barcodes)
    present=barcodes.isin(table.BARCODE)
    eligible=present & aligned['DROPLET.TYPE'].eq('SNG').to_numpy() & aligned['SNG.POSTERIOR'].gt(min_singlet_probability).to_numpy()
    donors=aligned.loc[eligible,'SNG.BEST.GUESS']
    if donors.isna().any() or not set(donors).issubset(expected):
        raise ValueError('Eligible assignment has an unexpected donor identity')
    aligned['assignment_present']=present
    aligned['assignment_eligible']=eligible
    aligned['assignment_exclusion_reason']=np.select(
        [~present,aligned['DROPLET.TYPE'].ne('SNG'),aligned['SNG.POSTERIOR'].le(min_singlet_probability)],
        ['missing_assignment','not_source_singlet','singlet_probability_at_or_below_threshold'],default='')
    aligned.index.name='BARCODE'
    return aligned
