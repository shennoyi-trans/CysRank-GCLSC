"""Generate true single-Cys insertion candidates and structure preparation requests."""
import argparse
from src.paths import root_path
from src import design


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    prepare = sub.add_parser("prepare", help="Enumerate insertion candidates; no fabricated structure scores")
    prepare.add_argument("--sequence", required=True)
    prepare.add_argument("--record-id", default="peptide")
    prepare.add_argument("--output", required=True, type=root_path)
    score = sub.add_parser("score", help="Score insertion candidates using their full-atom PDB structures")
    score.add_argument("--candidates", required=True, type=root_path)
    score.add_argument("--structures", required=True, type=root_path)
    score.add_argument("--checkpoint", type=root_path, default=root_path("models/site_predictor.pt"))
    score.add_argument("--output", required=True, type=root_path)
    score.add_argument("--top-k", type=int, default=3)
    score.add_argument("--allow-partial", action="store_true")
    run = sub.add_parser("run", help="Sequence -> local ESMFold -> CysRank Top3 -> external handoff")
    run.add_argument("--sequence", required=True)
    run.add_argument("--record-id", default="peptide")
    run.add_argument("--output", required=True, type=root_path)
    run.add_argument("--checkpoint", type=root_path, default=root_path("models/site_predictor.pt"))
    run.add_argument("--top-k", type=int, default=3)
    run.add_argument("--fold-model", default="facebook/esmfold_v1", help="HF model ID or local model directory")
    run.add_argument("--device", choices=["auto", "cpu", "cuda"], default="auto")
    run.add_argument("--chunk-size", type=int, default=32)
    run.add_argument("--local-files-only", action="store_true")
    args = parser.parse_args()
    if args.command == "run":
        from src.folding import run_design
        run_design(args)
        return
    if args.command == "score":
        design.score(args)
        return
    rows = design.prepare(args.sequence, args.record_id, args.output)
    print(f"Prepared {len(rows)} insertion sites in {args.output}; candidate structures required before scoring.")


if __name__ == "__main__":
    main()
