"""Export simulation requests, train a separate barrier surrogate, or predict barriers."""
import argparse
from src.paths import root_path
from src import barrier


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subs = parser.add_subparsers(dest="command", required=True)
    request = subs.add_parser("requests", help="Export pending JSONL for the calculation team")
    request.add_argument("--input", required=True, type=root_path)
    request.add_argument("--output", required=True, type=root_path)
    request.set_defaults(run=barrier.requests)
    train = subs.add_parser("train", help="Train using reviewed external simulation results")
    train.add_argument("--input", required=True, type=root_path)
    train.add_argument("--output", required=True, type=root_path)
    train.add_argument("--alpha", type=float, default=1.0)
    train.add_argument("--pretrained", type=root_path, default=root_path("models/pretrained/v_48_020.pt"))
    train.set_defaults(run=barrier.train)
    predict = subs.add_parser("predict", help="Predict simulation barriers, not experimental outcomes")
    predict.add_argument("--model", required=True, type=root_path)
    predict.add_argument("--input", required=True, type=root_path)
    predict.add_argument("--protocol", required=True, type=root_path)
    predict.add_argument("--output", required=True, type=root_path)
    predict.add_argument("--pretrained", type=root_path, default=root_path("models/pretrained/v_48_020.pt"))
    predict.set_defaults(run=barrier.predict)
    args = parser.parse_args()
    args.run(args)


if __name__ == "__main__":
    main()
