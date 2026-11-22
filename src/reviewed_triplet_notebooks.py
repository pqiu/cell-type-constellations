"""Prepare explicitly reviewed specimen triplets without guessing clinical identities."""
import json
from scripts.download_pbmc_asset import ROOT, download
from scripts.prepare_pbmc_notebooks import notebook
from src.original_triplet_counts import audit


def source_member_url(spec, role, endings):
    """Resolve a published member, including an explicitly reviewed basename."""
    urls = [u.replace('ftp://', 'https://') for k, vs in spec['source_metadata'].items()
            if k.startswith('supplementary_file') for u in vs]
    overrides = spec.get('source_member_overrides', {})
    if set(overrides) - {'genes_member', 'barcodes_member', 'matrix_member'}:
        raise ValueError('Unknown original source member role')
    if role in overrides:
        name = overrides[role]
        if not name.startswith(spec['accession'] + '_') or '/' in name or '\\' in name:
            raise ValueError('Source member override must be an exact original basename')
        hits = [u for u in urls if u.rsplit('/', 1)[-1] == name]
    else:
        hits = [u for u in urls if u.endswith(endings)]
    if len(hits) != 1:
        raise ValueError('Ambiguous original RNA triplet')
    return hits[0]


def prepare(accession, specimens, title, evidence, publications, early_min_rna_qc=False, dataset_prefix=None):
    if not specimens or len({s['accession'] for s in specimens}) != len(specimens):
        raise ValueError('Nonempty unique reviewed GEO records required')
    groups = {}; seen = set(); audits = {}
    for spec in specimens:
        row = spec['source_metadata']; gsm = spec['accession']
        if row['geo_accession'] != [gsm] or row['organism_ch1'] != ['Homo sapiens']:
            raise ValueError('Original source record/organism mismatch')
        urls = [u.replace('ftp://', 'https://') for k, vs in row.items() if k.startswith('supplementary_file') for u in vs]
        paths = {}; spec['assets'] = []
        for role, endings in [('genes_member', ('features.tsv.gz', 'genes.tsv.gz')), ('barcodes_member', ('barcodes.tsv.gz',)), ('matrix_member', ('matrix.mtx.gz',))]:
            url = source_member_url(spec, role, endings); name = url.rsplit('/', 1)[1]
            if not name.startswith(gsm + '_'): raise ValueError('Source asset mismatch')
            paths[role] = download(url, ROOT / f'data/downloads/{accession}_assets' / name)
            spec[role] = name; spec['assets'].append(dict(url=url, filename=name))
        result = audit(paths, seen, allowed_non_rna_types=spec.get('allowed_non_rna_feature_types', ()),
            allow_plaintext_gz_members=spec.get('allow_plaintext_gz_members', False),
            exclude_ambiguous_barcodes=spec.get('exclude_ambiguous_barcodes', False),
            non_rna_feature_ids=spec.get('non_rna_feature_ids', ()))
        spec['source_count_review'] = result; audits[gsm] = result
        groups.setdefault(result['gene_coordinate_sha256'], []).append(spec)
        print(dict(specimen=gsm, donor=spec['donor_id'], cells=result['cells']), flush=True)
    export_groups(accession, specimens, title, evidence, publications, groups, audits, early_min_rna_qc=early_min_rna_qc, dataset_prefix=dataset_prefix)


def export_groups(accession, specimens, title, evidence, publications, groups, audits, early_min_rna_qc=False, dataset_prefix=None):
    """Enable a notebook only after all original count and coordinate audits pass."""
    prefix = dataset_prefix or accession
    if prefix != accession:
        import re
        if not re.fullmatch(re.escape(accession) + r'_[a-z][a-z0-9_]*', prefix):
            raise ValueError('Compartment output prefix must retain original study accession')
    queue = []
    for number, (sha, samples) in enumerate(sorted(groups.items()), 1):
        cfg = json.loads((ROOT / 'references/datasets/configs/GSE136103.json').read_text())
        cfg.update(dataset_id=f'{prefix}_genes{number:03d}', study_id=accession, geo_samples=samples,
            sample_mapping_evidence=evidence,
            source_metadata=dict(accession=accession, publications=publications,
                source_url='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=' + accession,
                access='Public original RNA count matrices', license='not_specified_in_GEO',
                access_terms_url='https://www.ncbi.nlm.nih.gov/home/about/policies/',
                gene_coordinate_sha256=sha, source_processing=samples[0]['source_metadata']['data_processing']))
        if early_min_rna_qc:
            cfg['preconversion_min_rna_counts'] = cfg['qc']['min_counts']
            cfg['source_metadata']['source_barcode_scope'] = 'Original raw droplet whitelist; barcode rows are not called cells. Existing minimum RNA QC applied before constructing AnnData; no new threshold or RNA count changes.'
        path = ROOT / f'references/datasets/configs/{cfg["dataset_id"]}.json'
        if (ROOT / f'results/star_maps/{cfg["dataset_id"]}/run.json').exists(): raise ValueError('Refuse completed config changes')
        path.write_text(json.dumps(cfg, indent=2)); queue.append(path.relative_to(ROOT).as_posix())
    (ROOT / f'references/datasets/{prefix}_processing_queue.json').write_text(json.dumps(queue, indent=2))
    (ROOT / f'references/datasets/{prefix}_original_count_audits.json').write_text(json.dumps(audits, indent=2))
    nb = notebook(title, queue[0], True)
    if prefix == accession:
        nb['cells'][5]['source'] = ['from scripts.grouped_geo_notebooks import run\n', f'reports=run("{accession}")\n']
    else:
        nb['cells'][5]['source'] = ['import torch\n', 'from scimilarity import CellAnnotation\n',
            'from src.pbmc_pipeline import process\n', 'torch.set_num_threads(8); torch.manual_seed(0)\n',
            'annotation=CellAnnotation(str(ROOT / "models/model_v1.1"), use_gpu=torch.cuda.is_available())\n',
            f'queue=json.loads((ROOT / "references/datasets/{prefix}_processing_queue.json").read_text())\n',
            'reports=[process(ROOT / path, annotation=annotation) for path in queue]\n']
    nb['cells'][6]['source'] = ['import pandas as pd\n', 'display(pd.DataFrame(reports))\n']
    nb['cells'][7]['source'] = ['# Each export validates before source-count cleanup.\n']
    import nbformat
    nbformat.validate(nbformat.from_dict(nb))
    (ROOT / f'notebooks/datasets/{prefix}.ipynb').write_text(json.dumps(nb, indent=1))
    print(dict(source_observations=len(specimens), donors=len({s['donor_id'] for s in specimens}), RNA_gene_groups=len(groups)), flush=True)

