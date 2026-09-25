"""Exploratory redox-feature sensitivity check; does NOT construct reduced coordinates."""
import argparse
import copy
import csv
import json
import statistics
import sys
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from Bio.PDB import PDBParser
from Bio.Data.PDBData import protein_letters_3to1
from src import core, pipeline, structure


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default='results/2mi1_reduction_case_20260924')
    args = parser.parse_args()
    out = Path(args.output)
    if not out.is_absolute():
        out = ROOT / out
    if out.exists():
        raise FileExistsError('Use a new output directory')
    pdb = ROOT / 'data/structures/test.pdb'
    checkpoint = ROOT / 'models/site_predictor.pt'
    hashes = {'pdb': core.sha(pdb), 'checkpoint': core.sha(checkpoint)}
    models = list(PDBParser(QUIET=True).get_structure('2mi1', str(pdb)))
    inputs, evidence = [], []
    for index, model in enumerate(models, 1):
        sequence = ''.join(protein_letters_3to1[r.resname] for r in model['A'])
        assert sequence == 'AGCKNFFWKTFTSC'
        merged = None
        for position in (3, 14):
            sample, details, _ = structure.features(pdb, sequence, position, 'A', index)
            assert details['target_pdb_residue']['number'] == position
            assert sample['is_disulfide'][position - 1] == 1
            evidence.append({'position_1based': position, **details})
            if merged is None:
                merged = sample
            else:
                assert sample['xyz'] == merged['xyz']
                for key in (*core.BIO, 'feature_mask'):
                    merged[key][position - 1] = sample[key][position - 1]
        merged['record_id'] = f'model_{index:02d}_oxidized'
        inputs.append(merged)
        reduced = copy.deepcopy(merged)
        reduced['record_id'] = f'model_{index:02d}_reduction_assumption'
        reduced['is_disulfide'][2] = reduced['is_disulfide'][13] = 0
        reduced['state_basis'] = 'User confirmed pre-reaction reduction on 2026-09-24; coordinates remain oxidized'
        for key in ('seq', 'xyz', 'sasa', 'prev_is_pro', 'nearby_atoms', 'feature_mask'):
            assert reduced[key] == merged[key]
        inputs.append(reduced)
    out.mkdir(parents=True)
    core.write_jsonl(out / 'model_inputs.jsonl', inputs)
    core.write_jsonl(out / 'structure_evidence.jsonl', evidence)
    pipeline.predict(SimpleNamespace(input=out/'model_inputs.jsonl', checkpoint=checkpoint, output=out/'scores.jsonl'))
    scores = {(r['record_id'], r['position_1based']): r for r in core.load_jsonl(out/'scores.jsonl')}
    comparisons = []
    for row in inputs:
        a, b = [scores[(row['record_id'], p)] for p in (3, 14)]
        comparisons.append({'model_index': int(row['record_id'].split('_')[1]),
            'scenario': 'oxidized' if row['record_id'].endswith('_oxidized') else 'reduction_assumption',
            'cys3_score': a['raw_score'], 'cys14_score': b['raw_score'],
            'cys3_minus_cys14': a['raw_score']-b['raw_score'],
            'cys3_higher': a['raw_score'] > b['raw_score'],
            'cys3_filter_pass': a['passes_hard_filter'], 'cys14_filter_pass': b['passes_hard_filter'],
            'cys3_classified_positive': a['predicted_positive'], 'cys14_classified_positive': b['predicted_positive'],
            'cys3_sasa_A2': row['sasa'][2], 'cys14_sasa_A2': row['sasa'][13],
            'cys3_nearby_atoms': row['nearby_atoms'][2], 'cys14_nearby_atoms': row['nearby_atoms'][13]})
    with (out/'comparison.csv').open('w', encoding='utf-8-sig', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(comparisons[0]))
        writer.writeheader(); writer.writerows(comparisons)
    summaries = {}
    for scenario in ('oxidized', 'reduction_assumption'):
        items = [r for r in comparisons if r['scenario'] == scenario]
        summaries[scenario] = {'conformers': len(items), 'cys3_higher_count': sum(r['cys3_higher'] for r in items)}
        for key in ('cys3_score', 'cys14_score', 'cys3_minus_cys14'):
            values = [r[key] for r in items]
            summaries[scenario][key] = {'mean': statistics.mean(values), 'min': min(values), 'max': max(values)}
    overlap = [r['seq'] for r in core.load_jsonl(ROOT/'data/raw/训练正例.jsonl') if r['seq'] in 'AGCKNFFWKTFTSC']
    manifest = {'model_sha256': hashes['checkpoint'], 'pdb_sha256': hashes['pdb'],
        'script_sha256': core.sha(__file__), 'structure_code_sha256': core.sha(structure.__file__),
        'model_version': core.torch.load(checkpoint, map_location='cpu', weights_only=True)['config']['version'],
        'scope': 'Existing Cys scoring, all PDB conformers; no insertion or model tuning',
        'experimental_statement': 'User reports pre-reaction reduction; oxidized state cannot react; expected Cys3 > Cys14 was known before scoring',
        'reduction_operation': 'Set is_disulfide=0 at Cys3 and Cys14 only; retain all oxidized coordinates and other features',
        'is_reduced_structure': False, 'independent_test': False, 'training_overlap_sequences': overlap,
        'geometry_failures': sum(not e['geometry_passed'] for e in evidence),
        'summaries': summaries,
        'limitations': ['Oxidized NMR conformers are not reduced-state conformers', 'No hydrogens added, geometry relaxation or reduced-state structure prediction performed',
            'SASA and crowding recomputed on full oxidized peptide with repository definitions; historical training feature definitions not independently verified',
            'Conformer variation is not a confidence interval or independent experimental replication', 'Classification score is not reaction efficiency; expected ordering known and source overlaps training']}
    core.write_json(out/'summary.json', manifest)
    lines = ['# 2MI1 已有 Cys 位点探索性评分', '',
        '实验条件更新：用户确认反应前做过还原，氧化态不能反应。本次不插入残基、不重训或调参。', '',
        '**本次是二硫键特征的敏感性检查，不是还原态结构预测。** 原 PDB 坐标保持不变；对照保留几何识别的 is_disulfide=1，实验还原假设组仅将两个 Cys 的该特征设为 0。SASA/邻近原子数均从完整氧化态肽重算，不能冒充还原态实测特征。', '',
        '## 全部构象的原始模型分数', '',
        '| 构象 | 氧化态 Cys3 | 氧化态 Cys14 | 还原假设 Cys3 | 还原假设 Cys14 | 还原假设排序 |',
        '| --- | --- | --- | --- | --- | --- |']
    for i in range(len(models)):
        a,b=comparisons[2*i:2*i+2]
        lines.append(f"| {i+1} | {a['cys3_score']:.6f} | {a['cys14_score']:.6f} | {b['cys3_score']:.6f} | {b['cys14_score']:.6f} | {'Cys3 > Cys14' if b['cys3_higher'] else 'Cys3 <= Cys14'} |")
    s = summaries['reduction_assumption']
    lines += ['', f"还原假设组：{s['cys3_higher_count']}/{s['conformers']} 个构象中 Cys3 分数更高；平均分 Cys3={s['cys3_score']['mean']:.6f}，Cys14={s['cys14_score']['mean']:.6f}。均值仅作构象汇总，不代表溶液态热力学加权。", '',
        '氧化态组两个位点都被二硫键硬规则排除；原始分数仍保留以便对照。还原假设组的规则通过情况与 0.5 分类阈值结果分别记录在 comparison.csv，不能把低分描述成结构排除。', '',
        '## 解释边界', '',
        f"当前训练含 {', '.join(overlap)}，与测试肽重叠；本结果不是独立验证，也不是两个位点效率的定量预测。实验排序在运行前已知，不能据此调参或筛选构象。",
        '结构状态仍不匹配：改变一个标记没有松开原二硫键约束形成的构象。即使排序一致，也只构成该近似下的一致性观察；不一致亦不能单独区分模型问题与还原态结构近似问题。', '',
        '后续更严格检验需要状态匹配的还原态结构，以及排除该母来源组的模型。当前正式模型和原始 PDB 未改动。', '',
        f"模型 SHA-256：`{hashes['checkpoint']}`。",
        '输入、40 个位点评分、特征证据和统计见同目录 JSONL/CSV/JSON。', '',
        '复现命令：`python tools/score_2mi1_reduction_case.py --output results/2mi1_reduction_case_new_run`。']
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    assert core.sha(pdb)==hashes['pdb'] and core.sha(checkpoint)==hashes['checkpoint']
    print(json.dumps(summaries, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
