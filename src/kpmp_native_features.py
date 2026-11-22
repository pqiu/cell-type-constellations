"""Select explicitly labeled RNA rows from reviewed KPMP count products."""
import numpy as np
import pandas as pd


def select_rna_features(features):
    if features.shape[1] not in [2,3,6] or not features[0].is_unique:
        raise ValueError('Unreviewed feature table schema or repeated source IDs')
    if features.shape[1]==2:
        if features[0].str.match(r'^chr[^:]+:[0-9]+').any():
            raise ValueError('Genomic intervals cannot be treated as unlabeled RNA genes')
        keep=np.ones(len(features),dtype=bool)
    else:
        if not set(features[2]).issubset({'Gene Expression','Peaks'}):
            raise ValueError('Unknown modality in source feature_type')
        keep=features[2].eq('Gene Expression').to_numpy()
        if not keep.any():raise ValueError('No explicitly labeled RNA features')
    chosen=features.loc[keep]
    var=pd.DataFrame({'feature_name':chosen[1].to_numpy()},index=chosen[0].to_numpy())
    if features.shape[1]>=3:var['source_feature_type']=chosen[2].to_numpy()
    if features.shape[1]==6:
        for column,name in [(3,'source_chromosome'),(4,'source_start'),(5,'source_end')]:var[name]=chosen[column].to_numpy()
    return keep,var
