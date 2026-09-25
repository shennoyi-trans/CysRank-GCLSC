"""Read-only model audit; diagnostic refits are never saved as production weights."""
import json
import statistics
from collections import Counter, defaultdict
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import torch
import numpy as np
from src import core,pipeline

def main():
    torch.set_num_threads(4)
    out=ROOT/'results/training_audit_20260924';out.mkdir(exist_ok=False)
    saved=torch.load(ROOT/'models/site_predictor.pt',map_location='cpu',weights_only=True)
    config=saved['config'];device=torch.device(config['encoder_device'])
    backbone=core.make_backbone(ROOT/'models/pretrained/v_48_020.pt',device)
    assert all(torch.equal(v.cpu(),saved['backbone_state_dict'][k]) for k,v in backbone.state_dict().items())
    rows=pipeline.prepare(ROOT/'data/raw',ROOT/'data/provenance/parent_sequences.jsonl')
    assert rows==core.load_jsonl(ROOT/'data/processed/training.jsonl')
    assert all(core.sha(ROOT/'data/raw'/k)==v for k,v in config['source_sha256'].items())
    assert all(core.sha(ROOT/'src'/k)==v for k,v in config['code_sha256'].items())
    train=[];xs=[];ys=[]
    for r in rows:
        x=core.features(backbone,r,device)
        for i,m in enumerate(r['label_mask']):
            if m:
                xs.append(x[i]);ys.append(r['ca_label'][i])
                train.append({'id':r['record_id'],'seq':r['seq'],'position':i+1,'parent':r['parent_group'],
                              'parent_position':r['parent_start_1based']+i,'y':r['ca_label'][i],
                              **{k:r[k][i] for k in core.BIO}})
    x=torch.stack(xs);y=torch.tensor(ys,dtype=torch.float32)
    head=torch.nn.Linear(363,1);head.load_state_dict(saved['head_state_dict']);head.eval()
    scores=core.score_head(head,x,saved['columns'])
    old=core.load_jsonl(ROOT/'logs/retrain_20260924/training_predictions.jsonl')
    assert [(s['seq'],s['position'],s['y']) for s in train]==[(s['seq'],s['position_1based'],s['y']) for s in old]
    replay,_=core.fit(x,y,torch.ones(len(y)),config['epochs'],saved['columns'],config['seed'])
    groups=defaultdict(list)
    for r in train:groups[(r['parent'],r['parent_position'])].append(r)
    aliases=[rs for rs in groups.values() if len(rs)>1]
    case=ROOT/'results/2mi1_reduced_capped_ph75_20260924/analysis'
    inputs=core.load_jsonl(case/'model_inputs.jsonl')
    reference={(r['record_id'],r['position_1based']):r['raw_score'] for r in core.load_jsonl(case/'scores.jsonl')}
    cpu=backbone.cpu()
    segments={'structure':(0,128),'sequence':(128,359),'is_disulfide':(359,360),'sasa':(360,361),'prev_is_pro':(361,362),'nearby_atoms':(362,363)}
    w=head.weight.detach()[0];b=float(head.bias.detach()[0]);details=[];case_x=[]
    for r in inputs:
        if '_production_' not in r['record_id']:continue
        f=core.features(cpu,r,torch.device('cpu'))
        for pos in [3,14]:
            xi=f[pos-1];case_x.append(xi)
            contribution={k:float((xi[lo:hi]*w[lo:hi]).sum()) for k,(lo,hi) in segments.items()}
            score=float(torch.sigmoid(xi@w+b))
            assert abs(score-reference[r['record_id'],pos])<2e-6
            details.append({'id':r['record_id'],'position':pos,'score':score,'bias':b,**contribution,
                            'bio':{k:r[k][pos-1] for k in core.BIO}})
    delta={k:statistics.mean(details[i+1][k]-details[i][k] for i in range(0,len(details),2)) for k in segments}
    sequence_deltas=[]
    xx=torch.stack(case_x)
    for offset in range(-5,6):
        lo=128+(offset+5)*21
        value=float(((xx[1::2,lo:lo+21]-xx[0::2,lo:lo+21])*w[lo:lo+21]).sum(1).mean())
        sequence_deltas.append({'offset':offset,'cys14_minus_cys3_logit':value})
    diagnostics={}
    variants={'same_training_replay':list(range(363)), 'without_biochemical':list(range(359)),
              'sequence_only':list(range(128,359)), 'structure_only':list(range(128)), 'biochemical_only':list(range(359,363))}
    for name,cols in variants.items():
        h,_=core.fit(x,y,torch.ones(len(y)),50,cols,config['seed'])
        p=core.score_head(h,xx,cols)
        diagnostics[name]={'train_accuracy':core.metrics(ys,core.score_head(h,x,cols))['accuracy'],
                           'cys3_higher':int((p[::2]>p[1::2]).sum()),'cys3_mean':float(p[::2].mean()),'cys14_mean':float(p[1::2].mean())}
    keep=[i for i,r in enumerate(train) if r['parent']!='old3_positive_05']
    h,_=core.fit(x[keep],y[keep],torch.ones(len(keep)),50,list(range(363)),config['seed'])
    p=core.score_head(h,xx,list(range(363)))
    diagnostics['exclude_test_parent']={'train_n':len(keep),'cys3_higher':int((p[::2]>p[1::2]).sum()),'cys3_mean':float(p[::2].mean()),'cys14_mean':float(p[1::2].mean())}
    stats={}
    for label in [0,1]:
        rr=[r for r in train if r['y']==label]
        stats[str(label)]={'n':len(rr),'parents':dict(Counter(r['parent'] for r in rr)),
                         'disulfide':dict(Counter(r['is_disulfide'] for r in rr)),
                         'lengths':[len(r['seq']) for r in rr],
                         'feature_ranges':{k:[min(r[k] for r in rr),max(r[k] for r in rr)] for k in core.BIO}}
    result={'source_code_data_hashes_match_training':True,'prepared_data_identical':True,'frozen_backbone_matches':True,
            'training_score_max_error':float(max(abs(float(s)-r['raw_score']) for s,r in zip(scores,old))),
            'retrain_weight_max_error':float((replay.weight-head.weight).abs().max()),
            'model_sha256':core.sha(ROOT/'models/site_predictor.pt'),'training_sites':train,'class_summary':stats,
            'same_parent_site_multiple_records':aliases,'bias':b,'bio_weights_normalized':dict(zip(core.BIO,w[359:].tolist())),
            'mean_cys14_minus_cys3_logit_contributions':delta,'sequence_offset_contributions':sequence_deltas,
            'diagnostic_refits_not_model_selection':diagnostics}
    core.write_json(out/'audit.json',result);core.write_jsonl(out/'case_attribution.jsonl',details)
    print(json.dumps({k:v for k,v in result.items() if k not in ['training_sites','same_parent_site_multiple_records']},ensure_ascii=False,indent=2))

if __name__=='__main__':main()
