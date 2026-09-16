"""One-command training, optimizer updates, inference and top-k screening demo."""
import argparse
import json
import time
from types import SimpleNamespace
from src.paths import root_path
from src.pipeline import train, core
from predict import export


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=root_path, default=root_path('logs/demo'))
    parser.add_argument('--epochs', type=int, default=50)
    parser.add_argument('--track', default='未指定（演示）')
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Choose a new output directory')
    config = json.loads(root_path('configs/train.json').read_text(encoding='utf-8'))
    for key in ('source', 'parents', 'pretrained'):
        config[key] = root_path(config[key])
    config.update(output=args.output / 'training', epochs=args.epochs, feedback=None)
    started = time.perf_counter()
    train(SimpleNamespace(**config))
    export(root_path('data/examples/input.jsonl'), args.output / 'training/final_model/site_predictor.pt',
           args.output / 'results.csv', args.track)
    core.write_json(args.output / 'timing.json', {'elapsed_seconds': time.perf_counter() - started,
                                                'epochs': args.epochs})


if __name__ == '__main__':
    main()
