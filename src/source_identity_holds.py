"""Keep numerically validated observations with unresolved identity out of eligibility."""
import hashlib
import json


def load_holds(root):
    path=root / 'references/datasets/source_identity_holds.json'
    if not path.exists():return []
    holds=json.loads(path.read_text(encoding='utf-8'))
    ids=[row['sample_id'] for row in holds]
    if len(ids)!=len(set(ids)):raise ValueError('Duplicate source identity hold identifiers')
    for row in holds:
        if row['status']!='unresolved' or not row.get('reason'):raise ValueError('Explicit unresolved hold and reason required')
        evidence=(root / row['evidence']).resolve()
        if not evidence.is_relative_to(root.resolve()):raise ValueError('Hold evidence must remain within project')
        if hashlib.sha256(evidence.read_bytes()).hexdigest()!=row['evidence_sha256']:
            raise ValueError('Source identity hold evidence changed')
    return holds


def annotate_holds(samples,holds):
    ids=[row['sample_id'] for row in holds]
    if len(ids)!=len(set(ids)):raise ValueError('Duplicate source identity hold identifiers')
    by_id={row['sample_id']:row for row in holds}
    result=samples.copy()
    result['source_identity_hold']=result.sample_id.isin(by_id)
    result['source_identity_hold_details']=result.sample_id.map(lambda sid:json.dumps(by_id[sid],sort_keys=True) if sid in by_id else '')
    return result
