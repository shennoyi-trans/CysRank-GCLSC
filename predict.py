"""Generate UTF-8 candidate CSV and the full site-scoring audit JSONL."""
import argparse
import csv
import hashlib
import json
import tempfile
from pathlib import Path
from types import SimpleNamespace
from src.paths import ROOT, root_path
from src.pipeline import predict, core

FIELDS = ['candidate_id', 'track', 'sequence', 'position_1based', 'prediction_score',
          'rank', 'model_version', 'model_sha256', 'run_version', 'structure_file',
          'is_disulfide', 'sasa_A2', 'prev_is_pro', 'nearby_atoms_5A', 'remarks']


def export(input_path, checkpoint, output, track, top_k=3):
    """Rank eligible existing sites globally; never invent inserted sequences."""
    if top_k < 1:
        raise ValueError('top-k must be positive')
    audit = output.with_suffix('.sites.jsonl')
    metadata = output.with_suffix('.run.json')
    for path in (output, audit, metadata):
        if path.exists():
            raise FileExistsError(path)
    rows = core.load_jsonl(input_path)
    by_id = {}
    for row in rows:
        record_id = row.get('record_id')
        if not isinstance(record_id, str) or not record_id or record_id in by_id:
            raise ValueError('Each input requires a unique nonempty record_id')
        # A structure link must refer to a real repository file when supplied.
        if row.get('structure_file'):
            structure = root_path(row['structure_file'])
            structure.relative_to(ROOT)
            if not structure.is_file():
                raise FileNotFoundError(structure)
        by_id[record_id] = row
    with tempfile.TemporaryDirectory() as directory:
        temporary = Path(directory) / 'sites.jsonl'
        predict(SimpleNamespace(input=input_path, checkpoint=checkpoint, output=temporary))
        scores = core.load_jsonl(temporary)
    # Fixed threshold from the delivered model; no tuning on evaluation data.
    eligible = [site for site in scores if site['predicted_positive']]
    eligible.sort(key=lambda site: (-site['raw_score'], site['record_id'], site['position_1based']))
    digest = hashlib.sha256(checkpoint.read_bytes()).hexdigest()
    model_version = core.torch.load(checkpoint, map_location='cpu', weights_only=True)['config']['version']
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        for rank, site in enumerate(eligible[:top_k], 1):
            row = by_id[site['record_id']]
            i = site['position_1based'] - 1
            structure = root_path(row['structure_file']).relative_to(ROOT).as_posix() if row.get('structure_file') else ''
            writer.writerow(dict(candidate_id=f"{site['record_id']}_C{i+1}", track=track,
                sequence=site['seq'], position_1based=i+1,
                prediction_score=f"{site['raw_score']:.8f}", rank=rank,
                model_version=model_version, model_sha256=digest, run_version='submission-v1',
                structure_file=structure, is_disulfide=row['is_disulfide'][i],
                sasa_A2=row['sasa'][i], prev_is_pro=row['prev_is_pro'][i],
                nearby_atoms_5A=row['nearby_atoms'][i],
                remarks='Existing Cys site; uncalibrated score, not reaction efficiency; '
                        + ('structure supplied by input provider' if structure else 'full-atom structure unavailable')
                        + ('; training-derived demonstration, not independent validation'
                           if input_path == root_path('data/examples/input.jsonl') else '')))
    core.write_jsonl(audit, scores)
    core.write_json(metadata, dict(run_version='submission-v1', model_version=model_version,
        model_sha256=digest, input_sha256=core.sha(input_path), track=track, top_k=top_k,
        threshold=0.5, scored_sites=len(scores), eligible_sites=len(eligible),
        exported_candidates=min(top_k, len(eligible)),
        sorting='score descending, record_id ascending, 1-based position ascending',
        uncertainty='Uncalibrated model scores; no validated confidence intervals',
        torch_version=str(core.torch.__version__), numpy_version=core.np.__version__))
    print(f'Exported {min(top_k, len(eligible))} candidates to {output.name}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=root_path, default=root_path('data/examples/input.jsonl'))
    parser.add_argument('--checkpoint', type=root_path, default=root_path('models/site_predictor.pt'))
    parser.add_argument('--output', type=root_path, default=root_path('results/results.csv'))
    parser.add_argument('--track', default='未指定（演示）')
    parser.add_argument('--top-k', type=int, default=3)
    args = parser.parse_args()
    export(args.input, args.checkpoint, args.output, args.track, args.top_k)


if __name__ == '__main__':
    main()
