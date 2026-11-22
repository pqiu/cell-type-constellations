"""Explicitly keep RNA in original 10x CITE-seq feature tables."""
from src.kpmp_native_features import select_rna_features

def select_cite_rna_features(features):
    if features.shape[1]!=3 or not set(features[2]).issubset({'Gene Expression','Antibody Capture'}):
        raise ValueError('Unreviewed original CITE feature modality/layout')
    keep=features[2].eq('Gene Expression').to_numpy()
    if not keep.any():raise ValueError('No original RNA features')
    selected,var=select_rna_features(features.loc[keep].reset_index(drop=True))
    assert selected.all()
    return keep,var
