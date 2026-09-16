"""Train the ZQY8 site head; all relative paths are repository-root-relative."""
import argparse
import json
from types import SimpleNamespace
from src.paths import root_path
from src.pipeline import train


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=root_path, default=root_path('configs/train.json'))
    parser.add_argument('--output', type=root_path)
    parser.add_argument('--epochs', type=int)
    parser.add_argument('--seed', type=int)
    parser.add_argument('--feedback', type=root_path)
    args = parser.parse_args()
    config = json.loads(args.config.read_text(encoding='utf-8'))
    for name in ('source', 'parents', 'pretrained', 'output'):
        config[name] = root_path(config[name])
    for name in ('output', 'epochs', 'seed'):
        if getattr(args, name) is not None:
            config[name] = getattr(args, name)
    config['feedback'] = args.feedback
    train(SimpleNamespace(**config))


if __name__ == '__main__':
    main()
