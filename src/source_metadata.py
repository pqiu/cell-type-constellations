"""Explicit source-field aliases; original metadata columns remain intact."""


def apply_metadata_mapping(obs, mapping):
    result = obs.copy()
    for target, source in mapping.get('aliases', {}).items():
        if target in result.columns:
            raise ValueError('Refuse to overwrite original metadata: ' + target)
        if source not in result.columns:
            raise ValueError('Missing source metadata field: ' + source)
        result[target] = result[source].copy()
    for target, value in mapping.get('constants', {}).items():
        if target in result.columns:
            existing = result[target]
            matches = existing.isna().all() if value is None else existing.eq(value).all()
            if not matches:
                raise ValueError('Refuse to overwrite original metadata: ' + target)
            # An identical source-native constant is already present.  Treat the
            # mapping as idempotent while preserving the original column.
            continue
        result[target] = value
    for target, rule in mapping.get('value_maps', {}).items():
        if target in result.columns:
            raise ValueError('Refuse to overwrite original metadata: ' + target)
        source = rule['source']
        if source not in result.columns:
            raise ValueError('Missing source metadata field: ' + source)
        values = result[source].astype(str).map(rule['values'])
        if values.isna().any():
            raise ValueError('Unmapped source identity: ' + source)
        result[target] = values
    return result
