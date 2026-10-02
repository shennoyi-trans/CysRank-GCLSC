"""Isolated raw-value fitting of draft feedback; no asserted barrier units/protocol.

Uses the existing frozen ProteinMPNN encoder and ridge implementation. Does not
change the production feedback validator, reviewed flags, or classifier weights.
"""
import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
from src import barrier, core, pipeline, structure


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--alpha', type=float, default=1.0)
    args = parser.parse_args()
    source = args.input.resolve()
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(output)
    pretrained = ROOT / 'models/pretrained/v_48_020.pt'
    classifier = ROOT / 'models/site_predictor.pt'
    unchanged = {str(p):core.sha(p) for p in (source,pretrained,classifier)}
    rows = core.load_jsonl(source)
    if len(rows) < 2:
        raise ValueError('At least two distinct draft observations required')
    identities, sites, target, validation = set(), set(), [], []
    for row in rows:
        if row['evidence_type'] != 'simulation' or row['outcome'] is not None or row['split'] != 'train':
            raise ValueError('Trial expects only simulation train observations with no experimental outcome')
        if row['observation_id'] in identities:
            raise ValueError('Duplicate observation')
        identities.add(row['observation_id'])
        key = (row['sample']['seq'],row['position_1based'])
        if key in sites:
            raise ValueError('Duplicate sequence/site')
        sites.add(key)
        core.validate(row['sample'])
        if row['barrier']['sample_sha256'] != barrier.sample_digest(row['sample']):
            raise ValueError('Sample hash mismatch')
        # Explicitly unknown units are retained. This trial must not mix known units.
        if row['barrier']['unit'] is not None or row['barrier']['quantity'] is not None:
            raise ValueError('This draft-only script expects unknown unit/quantity; use formal training for completed feedback')
        value = row['barrier']['value']
        barrier.finite(value,'raw target')
        original_path = Path(row['provenance'])
        if core.sha(original_path) != row['completion']['source_sha256']:
            raise ValueError('Original return file changed')
        matches = [r for r in core.load_jsonl(original_path) if r['seq_id']==row['completion']['source_seq_id']]
        if len(matches)!=1 or matches[0]['seq']!=row['sample']['seq']:
            raise ValueError('Source record mismatch')
        if matches[0]['energy_barriers'][row['position_1based']-1] != value:
            raise ValueError('Source target mismatch')
        pdb = Path(row['completion']['model_input_structure_reference'])
        measured,evidence,_ = structure.features(pdb,row['sample']['seq'],row['position_1based'],'A')
        if not evidence['geometry_passed']:
            raise ValueError('Unresolved structural geometry problem')
        if any(measured[k]!=row['sample'][k] for k in ('seq',*core.BIO,'xyz')):
            raise ValueError('PDB and model input disagree')
        try:
            pipeline.validate_event(row)
            formal_error = None
        except ValueError as exc:
            formal_error = str(exc)
        validation.append({'observation_id':row['observation_id'],'geometry_passed':True,
                           'pdb_sha256':core.sha(pdb),'source_value_preserved':True,
                           'formal_validation_error':formal_error,
                           'formal_training_eligible':barrier.eligible(row) and formal_error is None})
        target.append(value)
    print(f'Extracting existing 363-dimensional features for {len(rows)} observations',flush=True)
    x = barrier.vectors(rows,pretrained)
    y = np.array(target,dtype=np.float64)
    if x.shape != (len(rows),363) or not np.isfinite(x).all():
        raise ValueError('Invalid feature matrix')
    head = barrier.fit_ridge(x,y,args.alpha)
    pred = barrier.apply_regressor(head,x)
    if not np.isfinite(pred).all():
        raise ValueError('Nonfinite predictions')
    error = pred-y
    model = {'version':'draft-feedback-raw-value-trial-v1','kind':'frozen_ProteinMPNN_ridge_raw_value_trial',
             'production_ready':False,'unit':None,'quantity':None,'protocol_confirmed':False,
             'scope':'Experimental numeric fit to supplied simulation values in unknown source units; not a validated barrier predictor',
             'alpha':args.alpha,'head':head,'input_sha256':core.sha(source),
             'pretrained_sha256':core.sha(pretrained),'training_observations':[r['observation_id'] for r in rows],
             'training_parent_groups':sorted({r['parent_group'] for r in rows}),
             'code_sha256':{'runner':core.sha(__file__),'barrier':core.sha(barrier.__file__),'core':core.sha(core.__file__)}}
    predictions = [{'seq':r['sample']['seq'],'position_1based':r['position_1based'],
                    'observation_id':r['observation_id'],'evaluation':'training_set_replay',
                    'target_raw_value':float(t),'predicted_raw_value':float(p),
                    'absolute_error_raw_value':float(abs(p-t)),'unit':None}
                   for r,t,p in zip(rows,y,pred)]
    summary = {'status':'completed_exploratory_fit','training_observations':len(rows),
               'parent_groups':len(model['training_parent_groups']),'feature_dimension':x.shape[1],
               'alpha':args.alpha,'target_unit':None,'target_quantity':None,
               'training_mae_raw_value':float(np.mean(np.abs(error))),
               'training_rmse_raw_value':float(np.sqrt(np.mean(error**2))),
               'training_mean_predictor_mae_raw_value':float(np.mean(np.abs(y-y.mean()))),
               'holdout':None,'validation_note':'All observations used to fit; training errors do not measure generalization',
               'formal_feedback_training_completed':False,'formal_model_ready':False,
               'draft_reason':'Unknown unit, quantity, laboratory protocol and review; original draft flags retained',
               'structure_note':'Two unrelaxed predictions and one locally minimized free-terminus model input; not laboratory calculation structures',
               'production_model_modified':False,'observations':validation}
    output.mkdir(parents=True)
    core.write_json(output/'trial_model.json',model)
    core.write_jsonl(output/'input.draft.jsonl',rows)
    core.write_jsonl(output/'training_predictions.jsonl',predictions)
    np.save(output/'training_features.npy',x,allow_pickle=False)
    reloaded=json.loads((output/'trial_model.json').read_text(encoding='utf-8'))
    replay=barrier.apply_regressor(reloaded['head'],np.load(output/'training_features.npy',allow_pickle=False))
    np.testing.assert_allclose(replay,pred,rtol=0,atol=1e-12)
    if any(core.sha(Path(p))!=digest for p,digest in unchanged.items()):
        raise ValueError('Source or production weights changed during trial')
    summary['model_reload_verified']=True
    core.write_json(output/'summary.json',summary)
    (output/'executed_script.py').write_text(Path(__file__).read_text(encoding='utf-8'),encoding='utf-8')
    print(json.dumps({'summary':summary,'predictions':predictions,'output':str(output)},ensure_ascii=False,indent=2),flush=True)


if __name__=='__main__':
    main()
