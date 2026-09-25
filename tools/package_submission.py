"""Build a local review archive; never upload data. See the submission checklist."""
import argparse
import hashlib
import json
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

ROOT = Path(__file__).resolve().parents[1]
FILES = ('README.md', 'requirements.txt', 'requirements-md.txt', 'requirements-propka.txt', 'train.py', 'predict.py', 'screen.py',
         'design.py', 'feedback.py', 'barrier.py', 'run.py', 'app.py', 'requirements-folding.txt')
DIRECTORIES = ('configs', 'data', 'src', 'models', 'notebooks', 'results',
               'logs/retrain_20260924', 'logs/verification_20260924',
               'logs/retrain_sst_20260924',
               'logs/reduced_md_20260924', 'logs/organization_20260924',
               'docs', 'licenses', 'tests', 'tools', 'web')
EXCLUDED = {'__pycache__', '.ipynb_checkpoints', '.git', '.venv', 'tmp', 'folding-cache', 'web_jobs'}


def submission_files():
    files = [ROOT / name for name in FILES]
    for name in DIRECTORIES:
        directory = ROOT / name
        if not directory.is_dir():
            raise FileNotFoundError(directory)
        files.extend(p for p in directory.rglob('*') if p.is_file())
    selected = []
    for path in sorted(set(files)):
        if any(part in EXCLUDED for part in path.relative_to(ROOT).parts) or path.suffix in {'.pyc', '.pyo'}:
            continue
        # Refuse linked files/environments rather than packaging their targets.
        if path.is_symlink() or any(p.is_symlink() for p in path.parents if p != ROOT):
            raise ValueError(f'Symlink not permitted: {path}')
        path.resolve().relative_to(ROOT)
        if not path.is_file():
            raise FileNotFoundError(path)
        selected.append(path)
    return selected


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default='tmp/agents/submission/CysRank-GCLSC.zip')
    args = parser.parse_args()
    output = Path(args.output).expanduser()
    if not output.is_absolute():
        output = ROOT / output
    output = output.resolve()
    files = submission_files()
    if output in files or output.exists():
        raise FileExistsError('Choose a new output archive outside submission inputs')
    output.parent.mkdir(parents=True, exist_ok=True)
    manifest = []
    # Exclusive creation protects previous packages.
    with ZipFile(output, 'x', compression=ZIP_DEFLATED) as archive:
        for path in files:
            content = path.read_bytes()
            relative = path.relative_to(ROOT).as_posix()
            archive.writestr('CysRank-GCLSC/' + relative, content)
            manifest.append({'path': relative, 'bytes': len(content),
                             'sha256': hashlib.sha256(content).hexdigest()})
        archive.writestr('CysRank-GCLSC/SUBMISSION_MANIFEST.json',
                         json.dumps(manifest, ensure_ascii=False, indent=2) + '\n')
    print(f'Packaged {len(manifest)} files: {output}')


if __name__ == '__main__':
    main()
