"""Fixed lambda=0.5 PROPKA sensitivity analysis of the capped reduced 2MI1 case.

Uses PROPKA 3.5.1 with an explicit primary-amide input typing adapter.
No training, lambda search, structural selection or original score overwrite.
"""
import argparse
import csv
import hashlib
import io
import json
import math
import logging
import statistics
from importlib.metadata import version
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LAMBDA = 0.5
REFERENCE_PKA = 8.0
PH = 7.5


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def correct(score, pka):
    if not 0 < score < 1 or not math.isfinite(pka):
        raise ValueError('Finite pKa and strictly interior original score required')
    logit = math.log(score) - math.log1p(-score)
    corrected = logit - LAMBDA * (pka - REFERENCE_PKA)
    sigmoid = (1 / (1 + math.exp(-corrected)) if corrected >= 0
               else math.exp(corrected) / (1 + math.exp(corrected)))
    return logit, corrected, sigmoid, 1 / (1 + 10 ** (pka - PH))


def molecule(path):
    from propka.lib import loadOptions
    from propka.parameters import Parameters
    from propka.input import read_parameter_file, read_pdb, protein_precheck
    from propka.molecular_container import MolecularContainer
    from propka.bonds import BondMaker
    assert version('propka') == '3.5.1', 'Adapter verified for PROPKA 3.5.1 only'
    options = loadOptions(['case.pdb'])
    params = read_parameter_file(options.parameters, Parameters())
    mol = MolecularContainer(params, options)
    mol.name = Path(path).stem
    # ACE must precede ALA as ATOM to avoid PROPKA's false N+ assignment.
    # Coordinates/names/elements and all other records are unchanged.
    text = '\n'.join('ATOM  ' + line[6:] if line.startswith('HETATM') and line[17:20] == 'ACE'
                     else line for line in Path(path).read_text().splitlines()) + '\n'
    mol.conformations, mol.conformation_names = read_pdb(io.StringIO(text), params, mol)
    mol.top_up_conformations()
    protein_precheck(mol.conformations, mol.conformation_names)
    mol.version.setup_bonding_and_protonation(mol)
    for name in mol.conformation_names:
        conf = mol.conformations[name]
        atoms = {(a.res_num, a.name): a for a in conf.atoms}
        n = atoms[15, 'N']
        assert n.res_name == 'NHE' and n.type == 'hetatm'
        assert n.get_bonded_heavy_atoms() == [atoms[14, 'C']]
        assert atoms[14, 'O'] in atoms[14, 'C'].bonded_atoms
        assert atoms[0, 'C'] in atoms[1, 'N'].bonded_atoms
        assert not any(a.terminal for a in conf.atoms)
        # PROPKA's ligand typer requires TWO heavy neighbors for N.am;
        # primary -CONH2 has ONE. Explicit chemistry-based correction to
        # built-in neutral NAM group, before extracting groups/protonation.
        n.sybyl_type = 'N.am'
        BondMaker().add_pi_electron_information(mol)
    mol.extract_groups()
    for name in mol.conformation_names:
        mol.conformations[name].sort_atoms()
    mol.find_covalently_coupled_groups()
    mol.calculate_pka()
    conf = mol.conformations['1A']
    groups = {(g.atom.res_num, g.type): g for g in conf.groups}
    assert groups[0, 'BBC'].charge == 0
    assert groups[1, 'BBN'].charge == 0
    amide = groups[15, 'NAM']
    assert amide.charge == 0 and not amide.titratable
    assert len(amide.atom.get_bonded_elements('H')) == 2
    assert not any(g.type in ('N+', 'C-', 'N31') for g in conf.groups)
    pka = {}
    for pos in (3, 14):
        g = groups[pos, 'CYS']
        assert not g.exclude_cys_from_results
        assert not any(a.element == 'S' for a in g.atom.bonded_atoms)
        assert math.isfinite(g.pka_value)
        pka[pos] = g.pka_value
    audit = [{'label': g.label, 'type': g.type, 'charge': g.charge,
              'titratable': g.titratable, 'pka': g.pka_value,
              'bonds': [[a.res_num, a.name] for a in g.atom.bonded_atoms]}
             for g in conf.groups]
    return pka, audit, text


def write_json(path, data):
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', default='results/2mi1_reduced_capped_ph75_20260924')
    parser.add_argument('--output', default='results/2mi1_propka_lambda05_20260924')
    args = parser.parse_args()
    folder, out = ROOT / args.input, ROOT / args.output
    out.mkdir(parents=True, exist_ok=False)
    (out / 'inputs').mkdir()
    logger = logging.getLogger('propka')
    handler = logging.FileHandler(out / 'propka.log', encoding='utf-8')
    handler.setLevel(logging.WARNING)
    logger.addHandler(handler)
    logger.propagate = False
    protocol = json.loads((folder / 'protocol.json').read_text())
    assert sha(ROOT / 'logs/retrain_20260924/final_model/site_predictor.pt') == protocol['checkpoint_sha256']
    scores = [json.loads(s) for s in (folder / 'analysis/scores.jsonl').read_text().splitlines()]
    scoremap = {(s['record_id'], s['position_1based']): s['raw_score'] for s in scores}
    snaps = [json.loads(s) for s in (folder / 'snapshots.jsonl').read_text().splitlines()]
    rows, audits = [], []
    for snap in snaps:
        path = folder / snap['pdb']
        assert sha(path) == snap['sha256']
        identity = f"start_{snap['start']:02d}_{snap['stage']}_{snap['time_ps']:03d}ps"
        pkas, audit, adapted = molecule(path)
        (out / 'inputs' / (identity + '.pdb')).write_text(adapted, encoding='utf-8')
        audits.append({'record_id': identity, 'groups': audit})
        for pos in (3, 14):
            score = scoremap[identity, pos]
            logit, corrected, value, alpha = correct(score, pkas[pos])
            rows.append({'record_id': identity, 'start': snap['start'], 'stage': snap['stage'],
                         'time_ps': snap['time_ps'], 'position': pos, 'original_score': score,
                         'pka': pkas[pos], 'thiolate_fraction_ph75': alpha,
                         'original_logit': logit, 'corrected_logit': corrected,
                         'corrected_score': value, 'lambda': LAMBDA,
                         'structure_file': path.relative_to(ROOT).as_posix(),
                         'structure_sha256': snap['sha256']})
    with (out / 'results.csv').open('w', encoding='utf-8-sig', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    write_json(out / 'typing_audit.json', audits)
    def summarize(selected):
        a = [r for r in selected if r['position'] == 3]
        b = [r for r in selected if r['position'] == 14]
        result = {'frames': len(a), 'original_cys3_higher': sum(x['original_score'] > y['original_score'] for x,y in zip(a,b)),
                  'corrected_cys3_higher': sum(x['corrected_score'] > y['corrected_score'] for x,y in zip(a,b)),
                  'pka_cys3_lower': sum(x['pka'] < y['pka'] for x,y in zip(a,b))}
        for pos, group in [(3,a),(14,b)]:
            result[f'cys{pos}'] = {key: {'mean': statistics.mean(r[key] for r in group),
                                       'min': min(r[key] for r in group), 'max': max(r[key] for r in group)}
                                  for key in ['original_score','corrected_score','pka','thiolate_fraction_ph75']}
        return result
    prod = [r for r in rows if r['stage'] == 'production']
    assert len(prod) == 200 and len(rows) == 240
    summary = {'lambda': LAMBDA, 'reference_pka': REFERENCE_PKA, 'pH': PH, 'propka_version': version('propka'),
               'formula': 'sigmoid(logit(original_score) - 0.5*(pKa-8.0))',
               'model_sha256': sha(ROOT/'logs/retrain_20260924/final_model/site_predictor.pt'), 'script_sha256': sha(__file__),
               'source_scores_sha256': sha(folder/'analysis/scores.jsonl'),
               'production': summarize(prod),
               'per_start': {str(i):summarize([r for r in prod if r['start']==i]) for i in range(1,11)},
               'typing': 'ACE ATOM; primary NHE explicitly N.am (neutral NAM); hydrogens regenerated by PROPKA; all caps and Cys audited',
               'typing_warning': 'Built-in NAM expects one H plus N (secondary amide); primary NH2 has two H plus N and emits a count warning per frame. All three interaction atoms are retained; neutral charge, two H and carbonyl connection asserted. See propka.log. Primary-amide transfer remains an approximation.',
               'limitations': ['Heuristic uncalibrated correction; no lambda tuning', 'Correlated snapshots and short neutral-thiol MD; not protonation equilibrium sampling',
                              'Related training fragments: not independent validation', 'Primary-amide typing adapter uses internal PROPKA 3.5.1 API; not unmodified CLI output',
                              'pH affects reported thiolate fraction; the specified score correction has no explicit pH term']}
    write_json(out/'summary.json', summary)
    handler.close()
    logger.removeHandler(handler)
    print(json.dumps(summary['production'], indent=2), flush=True)


if __name__ == '__main__':
    main()
