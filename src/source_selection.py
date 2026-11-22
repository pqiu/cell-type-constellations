"""Exclude explicitly reviewed source identities without altering observations."""
import numpy as np


def primary_data_mask(obs, override=None):
    """Select primary rows, allowing only a documented duplicate component."""
    if override is None:
        if 'is_primary_data' not in obs:
            return np.ones(len(obs), dtype=bool)
        return obs['is_primary_data'].eq(True).to_numpy()
    if not isinstance(override, dict):
        raise ValueError('Component primary override must be a mapping')
    required = {'column', 'values', 'primary_dataset_id', 'evidence'}
    if set(override) != required:
        raise ValueError('Component primary override has missing or unexpected fields')
    column = override['column']
    values = override['values']
    if not isinstance(column, str) or not column or column not in obs:
        raise ValueError('Component primary override column is missing')
    if (
        not isinstance(values, list)
        or not values
        or any(not isinstance(value, str) or not value for value in values)
        or len(values) != len(set(values))
    ):
        raise ValueError('Component primary override values must be unique nonempty strings')
    for field in ('primary_dataset_id', 'evidence'):
        if not isinstance(override[field], str) or not override[field].strip():
            raise ValueError('Component primary override requires ' + field)
    selected = obs[column].astype(str).isin(values).to_numpy()
    if not selected.any():
        raise ValueError('Component primary override selected no source rows')
    return selected


def source_exclusion_mask(obs, exclusions):
    excluded = np.zeros(len(obs), dtype=bool)
    for column, values in exclusions.items():
        if not isinstance(values, list) or any(not isinstance(v, str) for v in values):
            raise ValueError('Exclusion identities must be a list of strings')
        if column == '__index__':
            if column in obs or not obs.index.is_unique:
                raise ValueError('Exact source-cell exclusion requires an unambiguous native index')
            if len(values) != len(set(values)) or set(values) - set(obs.index.astype(str)):
                raise ValueError('Exact source-cell exclusions must be unique and present')
            excluded |= obs.index.astype(str).isin(values)
            continue
        if column not in obs:
            raise ValueError('Missing exclusion identity column: ' + column)
        excluded |= obs[column].astype(str).isin(values).to_numpy()
    return excluded
