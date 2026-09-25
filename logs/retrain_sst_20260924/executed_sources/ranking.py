"""Experimental within-peptide ordering; simulated conformers carry no outcome labels."""
import json
import torch
from . import core


def load(path, binary_rows):
    event = json.loads(path.read_text(encoding='utf-8'))
    if event['evidence_type'] != 'wet_experiment_ordering' or event['split'] != 'train':
        raise ValueError('Only experimental TRAIN ordering is allowed')
    if event['observation_weight'] != 1.0 or not event['provenance'] or not event['conditions']:
        raise ValueError('One unit-weight observation with provenance and conditions required')
    preferred, other = event['preferred_position'], event['other_position']
    if preferred == other or any(type(i) is not int or not 1 <= i <= len(event['sequence']) or event['sequence'][i-1] != 'C' for i in (preferred,other)):
        raise ValueError('Two distinct Cys positions required')
    groups = {r['parent_group'] for r in binary_rows if r['seq'] in event['sequence'] or event['sequence'] in r['seq']}
    if len(groups) != 1:
        raise ValueError('Reconcile ranking parent with existing sequence-related groups')
    seen = set()
    for item in event['samples']:
        row = item['sample']; core.validate(row)
        if row['seq'] != event['sequence'] or row['record_id'] in seen:
            raise ValueError('Sequence mismatch or repeated conformer')
        if any(k in row for k in ['ca_label','label_mask','outcome']):
            raise ValueError('Ordering conformers must not carry inferred binary labels')
        seen.add(row['record_id'])
        if core.sha(core.ROOT/item['structure_file']) != item['structure_sha256']:
            raise ValueError('Structure hash mismatch')
    if not seen:
        raise ValueError('No conformers supplied')
    event['parent_group'] = next(iter(groups))
    return event


def encode(event, backbone, device):
    pairs = [core.features(backbone, r['sample'], device) for r in event['samples']]
    return {'preferred':torch.stack([f[event['preferred_position']-1] for f in pairs]),
            'other':torch.stack([f[event['other_position']-1] for f in pairs])}


def stats(head, pairs, columns):
    a = core.score_head(head, pairs['preferred'], columns)
    b = core.score_head(head, pairs['other'], columns)
    return {'conformers':len(a), 'preferred_higher':int((a>b).sum()),
            'preferred_mean':float(a.mean()), 'other_mean':float(b.mean()),
            'scope':'ONE experimental ordering; conformers are correlated augmentation, not independent tests'}
