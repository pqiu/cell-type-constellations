"""Keep experimental observation views distinct from source biological specimens."""
import pandas as pd
from urllib.parse import quote


def read_samples_preserving_identity(source,cfg):
    """Keep nominal source keys intact while reading numeric support normally."""
    return pd.read_csv(source,dtype={column:str for column in cfg.get('physical_sample_columns',[])})


def attach_physical_identity(obs,cfg):
    columns=cfg.get('physical_sample_columns')
    if columns is None:return obs
    if not columns or not cfg.get('physical_sample_mapping_evidence'):
        raise ValueError('Physical specimen identity needs explicit source evidence')
    result=obs.copy()
    encoding=cfg.get('physical_sample_value_encoding')
    if encoding not in (None,'percent'):
        raise ValueError('Unsupported physical specimen value encoding: '+str(encoding))
    for column in columns:
        if column not in result or result[column].isna().any():
            raise ValueError('Missing physical specimen identity: '+column)
        values=result[column].astype(str)
        if values.str.lower().isin(['','unknown','nan','na']).any():
            raise ValueError('Unknown or ambiguous physical specimen identity: '+column)
        if encoding is None and values.str.contains('|',regex=False).any():
            raise ValueError('Unknown or ambiguous physical specimen identity: '+column)
    if 'physical_specimen_uid' in result:
        raise ValueError('Refuse to overwrite original physical specimen identity')
    identity=result[columns].astype(str)
    if encoding=='percent':
        # Reversible escaping keeps literal source separators distinct from the
        # separator used to combine multiple identity columns.  The original
        # source columns remain unchanged in outputs.
        identity=identity.map(lambda value: quote(value,safe=''))
    result['physical_specimen_uid']=cfg['study_id']+':'+identity.agg('|'.join,axis=1)
    check_physical_identity(result)
    return result


def normalize_physical_identity(samples):
    result=samples.copy()
    if 'physical_specimen_uid' not in result:result['physical_specimen_uid']=result.sample_id
    else:result['physical_specimen_uid']=result.physical_specimen_uid.fillna(result.sample_id)
    check_physical_identity(result)
    return result


def check_physical_identity(frame):
    if frame.physical_specimen_uid.isna().any() or frame.physical_specimen_uid.astype(str).str.lower().isin(['','unknown','na','nan']).any():
        raise ValueError('Unknown physical specimen UID')
    if frame.groupby('sample_id',observed=True).physical_specimen_uid.nunique().gt(1).any():
        raise ValueError('One observation view maps to multiple physical specimens')
    if frame.groupby('physical_specimen_uid',observed=True).donor_uid.nunique().gt(1).any():
        raise ValueError('A physical specimen maps to multiple donors')


def specimen_summary(samples):
    frame=normalize_physical_identity(samples)
    result=dict(specimens=int(frame.physical_specimen_uid.nunique()),observation_views=int(frame.sample_id.nunique()),
                donors=int(frame.donor_uid.nunique()),cells=int(frame.retained_cells.sum()))
    if 'supported_stars' in frame:
        # Several selection views cannot be combined into an artificial map.
        best=frame.groupby('physical_specimen_uid',observed=True).supported_stars.max()
        result.update(multi_star=int(best.ge(2).sum()),no_stars=int(best.eq(0).sum()))
    return result
