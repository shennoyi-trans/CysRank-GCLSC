"""Score all snapshots of the fixed capped reduced-peptide protocol, retaining failures."""
import argparse
import copy
import csv
import json
import statistics
import sys
from pathlib import Path
from types import SimpleNamespace

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.Data.PDBData import protein_letters_3to1
from src import core, pipeline


def extract(path, identity):
    chain=PDBParser(QUIET=True).get_structure(identity,str(path))[0]['A']
    residues=list(chain)
    assert [r.resname for r in residues][0]=='ACE' and residues[-1].resname=='NHE'
    peptide=residues[1:-1]
    assert [r.id[1] for r in peptide]==list(range(1,15))
    seq=''.join(protein_letters_3to1[r.resname] for r in peptide)
    assert seq=='AGCKNFFWKTFTSC'
    assert {'N','HN1','HN2'}=={a.name for a in residues[-1]}
    assert all('HG' in peptide[p-1] for p in (3,14))
    # Use all peptide + cap heavy atoms for SASA and crowding; exclude explicit water/ions.
    for r in residues:
        for atom in list(r):
            if atom.element in ('H','D'):r.detach_child(atom.id)
    atoms=list(chain.get_atoms());coords=np.array([a.coord for a in atoms],dtype=float)
    assert np.isfinite(coords).all()
    ShrakeRupley(probe_radius=1.4,n_points=960).compute(chain,level='A')
    row={'record_id':identity,'seq':seq,'xyz':[r['CA'].coord.astype(float).tolist() for r in peptide],
         **{k:[0]*14 for k in (*core.BIO,'feature_mask')}}
    sg_distance=float(np.linalg.norm(peptide[2]['SG'].coord-peptide[13]['SG'].coord))
    violations=[]
    for i,r in enumerate(peptide,1):
        for a,b,lo,hi in [('N','CA',1.2,1.7),('CA','C',1.2,1.7),('C','O',1.0,1.5)]:
            d=float(np.linalg.norm(r[a].coord-r[b].coord))
            if not lo<=d<=hi:violations.append({'residue':i,'bond':f'{a}-{b}','distance_A':d})
    for prev,nxt in zip(residues[:-1],residues[1:]):
        d=float(np.linalg.norm(prev['C'].coord-nxt['N'].coord))
        if not 1.1<=d<=1.6:violations.append({'residue':nxt.id[1],'bond':'previous_C-N','distance_A':d})
    for position in (3,14):
        r=peptide[position-1];idx=position-1;sg=r['SG'].coord
        row['sasa'][idx]=float(r['CB'].sasa+r['SG'].sasa)
        row['nearby_atoms'][idx]=sum(float(np.linalg.norm(a.coord-sg))<=5 and a is not r['SG'] for a in atoms)
        # Preserve the current model's geometric threshold, not a new learned rule.
        # If a reduced-topology frame still has close SGs, disclose the false bonding proxy.
        row['is_disulfide'][idx]=int(sg_distance<2.5)
        row['prev_is_pro'][idx]=int(idx>0 and seq[idx-1]=='P')
        row['feature_mask'][idx]=1
        for a,b,lo,hi in [('CA','CB',1.2,1.8),('CB','SG',1.5,2.1)]:
            d=float(np.linalg.norm(r[a].coord-r[b].coord))
            if not lo<=d<=hi:violations.append({'residue':position,'bond':f'{a}-{b}','distance_A':d})
    ids=[residues.index(a.parent) for a in atoms]
    clashes=[]
    for i in range(len(atoms)):
        for j in range(i+1,len(atoms)):
            if abs(ids[i]-ids[j])>1 and np.linalg.norm(coords[i]-coords[j])<1.2:
                clashes.append([atoms[i].parent.id[1],atoms[i].name,atoms[j].parent.id[1],atoms[j].name])
    evidence={'record_id':identity,'sg_distance_A':sg_distance,'geometry_pass':not violations and not clashes,
        'bond_violations':violations,'severe_overlaps':clashes,
        'sasa_scope':'CB+SG; full capped peptide heavy atoms; water/ions excluded; probe 1.4 A, 960 points',
        'nearby_scope':'All capped peptide heavy atoms within 5 A of SG excluding target SG only',
        'encoder_scope':'14 amino-acid CA coordinates; ACE/NHE omitted from sequence encoder but included in scalar features',
        'reduced_topology':True,'geometric_disulfide_proxy':sg_distance<2.5}
    return row,evidence


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',default='results/2mi1_reduced_capped_ph75_20260924')
    args=p.parse_args();folder=ROOT/args.input
    protocol=json.loads((folder/'protocol.json').read_text())
    assert protocol['status']=='complete' and protocol['preparation_pH']==7.5
    out=folder/'analysis';out.mkdir(exist_ok=False)
    checkpoint=ROOT/'models/site_predictor.pt'
    assert core.sha(checkpoint)==protocol['checkpoint_sha256']
    snaps=core.load_jsonl(folder/'snapshots.jsonl')
    assert len(snaps)==120 and sum(s['stage']=='production' for s in snaps)==100
    inputs=[];evidences=[];metadata={}
    for s in snaps:
        path=folder/s['pdb'];assert core.sha(path)==s['sha256']
        identity=f"start_{s['start']:02d}_{s['stage']}_{s['time_ps']:03d}ps"
        row,evidence=extract(path,identity)
        inputs.append(row);evidences.append(evidence);metadata[identity]=s
    core.write_jsonl(out/'model_inputs.jsonl',inputs)
    core.write_jsonl(out/'structure_evidence.jsonl',evidences)
    pipeline.predict(SimpleNamespace(input=out/'model_inputs.jsonl',checkpoint=checkpoint,output=out/'scores.jsonl'))
    scores={(r['record_id'],r['position_1based']):r for r in core.load_jsonl(out/'scores.jsonl')}
    comparisons=[]
    for row,e in zip(inputs,evidences):
        s=metadata[row['record_id']];a,b=[scores[(row['record_id'],i)] for i in (3,14)]
        comparisons.append({**{k:s[k] for k in ('start','stage','time_ps','pdb','ca_rmsd_from_prepared_A')},
            'cys3_score':a['raw_score'],'cys14_score':b['raw_score'],'difference':a['raw_score']-b['raw_score'],
            'cys3_higher':a['raw_score']>b['raw_score'],'geometry_pass':e['geometry_pass'],
            'sg_distance_A':e['sg_distance_A'],'both_pass_filters':a['passes_hard_filter'] and b['passes_hard_filter'],
            'cys3_sasa_A2':row['sasa'][2],'cys14_sasa_A2':row['sasa'][13],
            'cys3_nearby_atoms':row['nearby_atoms'][2],'cys14_nearby_atoms':row['nearby_atoms'][13]})
    with (out/'comparison.csv').open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(comparisons[0]));w.writeheader();w.writerows(comparisons)
    def stats(rows):
        result={'frames':len(rows),'cys3_higher':sum(r['cys3_higher'] for r in rows),
            'geometry_failures':sum(not r['geometry_pass'] for r in rows),'filter_failures':sum(not r['both_pass_filters'] for r in rows)}
        for key in ('cys3_score','cys14_score','difference','sg_distance_A','ca_rmsd_from_prepared_A'):
            vals=[r[key] for r in rows]
            result[key]={'mean':statistics.mean(vals),'min':min(vals),'max':max(vals)} if vals else None
        return result
    prod=[r for r in comparisons if r['stage']=='production']
    summary={'production':stats(prod),'minimized':stats([r for r in comparisons if r['stage']=='minimized']),
        'valid_production':stats([r for r in prod if r['geometry_pass'] and r['both_pass_filters']]),
        'per_start':{str(i):stats([r for r in prod if r['start']==i]) for i in range(1,11)},
        'checkpoint_sha256':core.sha(checkpoint),'analysis_script_sha256':core.sha(__file__),
        'limitations':['100 ps per start is exploratory and not convergence evidence',
            'Time-correlated snapshots are not independent experiments; conformer count is not statistical sample size',
            'TFA and reaction buffer not explicitly modeled; 0.15 M NaCl approximation',
            'Initial protonation chosen at pH 7.5, fixed neutral thiols; no pKa or reactive chemistry model',
            'Current model trained on related fragments, no independent validation or efficiency prediction',
            'Caps included in SASA/crowding but unsupported in sequence encoder; historical feature compatibility uncertain',
            'Compared with previous trial both caps and geometry change, cannot isolate a reduction-only effect']}
    core.write_json(out/'summary.json',summary)
    s=summary['production'];mins=summary['minimized']
    lines=['# 还原态 Ac–AGCKNFFWKTFTSC–NH₂：结构重建与探索性评分','',
        '用户确认：N 端乙酰化、C 端酰胺化、Cys 巯基不保护，三氟乙酸盐，pH 7.5、300 K。',
        '采用 Amber ff14SB、TIP3P 水、ACE/NHE 端基；NHE 是一级酰胺，不是 N-甲基酰胺。移除原二硫键与全部旧氢，重建两个中性 Cys–SH，最小化后在 NVIDIA OpenCL 上运行 NVT 模拟。',
        '每个原始构象固定 20 ps 平衡、100 ps 采样，2 fs 步长，每 10 ps 保存一帧；10 个起点共 1 ns 生产采样，另有 0.2 ns 平衡。未按评分选择构象。',
        '0.15 M NaCl 为溶剂假设；TFA 含量和缓冲体系未知，未显式模拟 TFA。pH 仅用于初始质子化分配，不是恒 pH 模拟。','',
        '## 结果','',
        f"生产快照中 Cys3 > Cys14：**{s['cys3_higher']}/{s['frames']}**。平均分：Cys3 **{s['cys3_score']['mean']:.6f}**，Cys14 **{s['cys14_score']['mean']:.6f}**。",
        f"最小化结构：{mins['cys3_higher']}/{mins['frames']} 个起点为 Cys3 更高。",
        f"生产快照 SG–SG 距离 {s['sg_distance_A']['min']:.2f}–{s['sg_distance_A']['max']:.2f} Å；几何检查失败 {s['geometry_failures']} 帧，现有过滤规则不通过 {s['filter_failures']} 帧。所有帧均保留，合格子集统计见 summary.json。",'',
        '| 起点 | Cys3 平均分 | Cys14 平均分 | Cys3 更高的帧数 |','| --- | --- | --- | --- |']
    for i,r in summary['per_start'].items():
        lines.append(f"| {i} | {r['cys3_score']['mean']:.6f} | {r['cys14_score']['mean']:.6f} | {r['cys3_higher']}/{r['frames']} |")
    lines+=['','## 解释边界','',
        '这是短时松弛后的计算构象集合，不是已验证或已收敛的真实还原态结构。100 帧有时间相关性，不能当作 100 次独立实验或给出独立样本置信区间；均值未进行热力学重加权。',
        '重新计算的 SASA 和拥挤度包含 ACE/NHE 的重原子；序列编码器仍仅接收 14 个标准残基。与旧的未封端、未松弛尝试相比，端基和构象同时改变，不宜只归因于还原。',
        '当前模型接触过该肽的两个训练片段，分数不是环化效率。本次观察无论是否符合预期，都不足以证明独立泛化；没有调整模型、阈值或模拟参数来迎合已知实验排序。',
        '正式模型和原始 PDB 保持不变。后续应使用来源组留出的模型，确认反应溶剂，并延长、多重复检验构象采样。','',
        '## 文件与复现','',
        '`../start_*/` 保存带正确端基的还原态 PDB、显式溶剂初始结构、system/integrator XML、最终状态与检查点、能量日志。`../protocol.json` 和 `../runs.json` 保存完整参数、版本与随机种子。',
        '`comparison.csv` 为全部结构排序，`scores.jsonl` 为全部位点分数；`model_inputs.jsonl` 和 `structure_evidence.jsonl` 保留特征与核验依据。',
        '复现：`python tools/rebuild_2mi1_reduced.py --output results/new_reduced_run --ph 7.5 --temperature 300`，完成后运行 `python tools/analyze_2mi1_reduced_md.py --input results/new_reduced_run`。安装 OpenMM 8.6.1、PDBFixer 1.12.0 和项目运行依赖；正式构建不调用 PDBFixer，保留该安装仅用于早期输入检查。',
        '早期目录 `2mi1_reduced_md_20260924`、`2mi1_reduced_md_20260924_v2` 为准备失败尝试；`2mi1_reduced_capped_md_20260924` 为被新 pH 信息取代的未完成尝试，均不纳入本报告。']
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    assert core.sha(checkpoint)==protocol['checkpoint_sha256']
    assert core.sha(ROOT/'data/structures/test.pdb')==protocol['pdb_sha256']
    print(json.dumps(summary['production'],indent=2),flush=True)


if __name__=='__main__':main()
