"""Use only explicit source hash calls; no donor assignment for negatives/doublets."""
import pandas as pd


def select(frame, tag_fractions, donor, specimen):
    if not {'CellName','HashID'}.issubset(frame.columns):raise ValueError('Missing source hash fields')
    if frame.CellName.isna().any() or not frame.CellName.is_unique or frame.CellName.str.strip().eq('').any():raise ValueError('Invalid source cell identities')
    allowed=set(tag_fractions)|{'Negative','Doublet'}
    if frame.HashID.isna().any() or not set(frame.HashID).issubset(allowed):raise ValueError('Unknown source hash assignment')
    if set(frame.HashID)&set(tag_fractions)!=set(tag_fractions):raise ValueError('Expected source fraction missing')
    kept=frame.loc[frame.HashID.isin(tag_fractions)].copy().set_index('CellName')
    kept['biological_donor']=donor;kept['biological_specimen']=specimen
    kept['source_fraction']=kept.HashID.map(tag_fractions)
    return kept,dict(source_calls=len(frame),selected_cells=len(kept),
        negative_calls=int(frame.HashID.eq('Negative').sum()),doublet_calls=int(frame.HashID.eq('Doublet').sum()),
        source_fraction_support=kept.source_fraction.value_counts().to_dict(),
        confidence='Categorical source calls; no probability scores deposited or inferred')
