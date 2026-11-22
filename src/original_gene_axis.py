"""Check native row labels against an original feature list, including R aliases.

Seurat's unique feature labels can append .1, .2, etc to repeated gene symbols.
Only aliases reconstructed from the complete ordered original feature list are
accepted. This does not rename rows, merge genes, or infer measurement from zeros.
"""
def r_unique_names(names):
    reserved = set(names)
    used = set()
    result = []
    for name in names:
        value = name
        if value in used:
            suffix = 1
            while f'{name}.{suffix}' in reserved or f'{name}.{suffix}' in used:
                suffix += 1
            value = f'{name}.{suffix}'
        used.add(value)
        result.append(value)
    return result


def review_original_axis(native_names, axis):
    original = set(axis['gene_ids']) | set(axis['gene_symbols'])
    unique = r_unique_names(axis['gene_symbols'])
    alias = {name: {'original_symbol': symbol, 'original_gene_id': gene_id}
             for name, symbol, gene_id in zip(unique, axis['gene_symbols'], axis['gene_ids'])
             if name != symbol and name not in original}
    missing = [name for name in native_names if name not in original and name not in alias]
    used = {name: alias[name] for name in native_names if name in alias}
    return missing, used
