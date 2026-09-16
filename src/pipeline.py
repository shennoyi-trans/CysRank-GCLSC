"""ZQY8 site predictor and provenance-aware feedback replay (not insertion design)."""
import argparse
import copy
import json
from collections import Counter
from pathlib import Path

import numpy as np
import torch
from . import core

ROOT = core.ROOT
RULES = {"exclude_disulfide": True, "exclude_prev_pro": True,
         "sasa_hard_filter": False, "nearby_hard_filter": False,
         "prev_pro_basis": "user_adopted_conservative_hypothesis_not_proven_impossibility"}


def reasons(row, i, rules):
    return (["disulfide"] if rules["exclude_disulfide"] and row["is_disulfide"][i] else []) + (
        ["previous_residue_proline"] if rules["exclude_prev_pro"] and row["prev_is_pro"][i] else [])


def positive_filename(source):
    matches = [name for name in ("训练正例.jsonl.txt", "训练正例.jsonl", "训练正例.jsonl.jsonl") if (source / name).is_file()]
    if len(matches) != 1:
        raise ValueError(f"Need exactly one main positive file, found: {matches}")
    return matches[0]


def prepare(source, parents_path):
    parents = core.load_jsonl(parents_path)
    rows = []
    # Deliberately do not open the low-precision group or validation files.
    for filename, outcome in ((positive_filename(source), 1), ("训练集负例.jsonl", 0)):
        for line, raw in enumerate(core.load_jsonl(source / filename), 1):
            core.validate(raw)
            mask = raw["ca_label"]
            if not any(mask) or any(v and raw["seq"][i] != "C" for i, v in enumerate(mask)):
                raise ValueError(f"Invalid tested mask: {filename}:{line}")
            matches = [(j + 1, p["seq"].find(raw["seq"])) for j, p in enumerate(parents)
                       if raw["seq"] in p["seq"]]
            if len(matches) != 1:
                raise ValueError(f"Unresolved parent: {filename}:{line}")
            parent, start = matches[0]
            r = copy.deepcopy(raw)
            r.update(record_id=f"{'positive' if outcome else 'negative'}_{line:03d}",
                     source_file=filename, source_line=line, source_ca_label=mask.copy(),
                     tested_mask=mask.copy(), label_mask=mask.copy(),
                     ca_label=[outcome if v else 0 for v in mask],
                     confirmed_negative_mask=[int(v and not outcome) for v in mask],
                     experiment_confirmed=True, evidence_type="wet_experiment",
                     confirmation_basis="user_report; raw assay records not independently audited",
                     supervision_weight=1.0, quality_group="main_only",
                     parent_group=f"old3_positive_{parent:02d}", parent_start_1based=start + 1,
                     parent_assignment_basis="exact_subsequence; homology not audited")
            rows.append(r)
    return rows


def validate_event(event):
    """Explicit full-structure sample + site-specific observation; never infer labels."""
    for key in ("observation_id", "evidence_type", "split", "parent_group", "assay_id",
                "conditions", "provenance", "sample", "position_1based", "outcome", "reviewed"):
        if key not in event:
            raise ValueError(f"Missing feedback field: {key}")
    for key in ("observation_id", "parent_group", "assay_id", "provenance"):
        if not isinstance(event[key], str) or not event[key].strip():
            raise ValueError(f"Empty feedback field: {key}")
    if not isinstance(event["conditions"], dict) or not event["conditions"]:
        raise ValueError("Nonempty reaction/simulation conditions required")
    if event["evidence_type"] not in ("wet_experiment", "simulation"):
        raise ValueError("Unknown evidence_type")
    if event["split"] not in ("train", "holdout") or type(event["reviewed"]) is not bool:
        raise ValueError("Invalid split/reviewed")
    core.validate(event["sample"])
    i = event["position_1based"]
    if type(i) is not int or not 1 <= i <= len(event["sample"]["seq"]) or event["sample"]["seq"][i-1] != "C":
        raise ValueError("Feedback target must be a valid Cys position")
    outcome = event["outcome"]
    if outcome is not None and (type(outcome) is not int or outcome not in (0, 1)):
        raise ValueError("outcome must be 0/1/null; efficiency does not imply binary success")
    if event["evidence_type"] == "wet_experiment" and event["reviewed"] and outcome is None:
        raise ValueError("Reviewed experimental feedback needs an explicit outcome")
    efficiency = event.get("cyclization_efficiency_fraction")
    if efficiency is not None and (type(efficiency) not in (int, float) or not np.isfinite(efficiency) or not 0 <= efficiency <= 1):
        raise ValueError("Efficiency must be a fraction in [0,1], not percent")


def ingest(args):
    events = core.load_jsonl(args.input)
    seen = set()
    for event in events:
        validate_event(event)
        if event["observation_id"] in seen:
            raise ValueError("Duplicate observation_id")
        seen.add(event["observation_id"])
    if args.output.exists():
        raise FileExistsError(args.output)
    # Archive every observation; only reviewed wet TRAIN records can reach the loss.
    core.write_jsonl(args.output, events)
    print(json.dumps({"archived": len(events), "eligible_for_training": sum(
        e["evidence_type"] == "wet_experiment" and e["reviewed"] and e["split"] == "train" for e in events)}))


def merge_feedback(rows, path):
    rows = copy.deepcopy(rows)
    if path is None:
        return rows, {"observations": 0, "accepted": 0, "simulation": 0, "holdout": 0, "unreviewed": 0}
    events = core.load_jsonl(path)
    seen = set()
    used = {(r["seq"], i + 1) for r in rows for i, m in enumerate(r["label_mask"]) if m}
    stats = Counter(observations=len(events), accepted=0, simulation=0, holdout=0, unreviewed=0)
    held_groups = set()
    for event in events:
        validate_event(event)
        if event["observation_id"] in seen:
            raise ValueError("Duplicate observation_id")
        seen.add(event["observation_id"])
        seq = event["sample"]["seq"]
        related_groups = {r["parent_group"] for r in rows if seq in r["seq"] or r["seq"] in seq}
        if related_groups and related_groups != {event["parent_group"]}:
            raise ValueError("Feedback parent conflicts with sequence containment; reconcile provenance")
        if event["split"] == "holdout":
            held_groups.add(event["parent_group"])
            stats["holdout"] += 1
            continue
        if event["evidence_type"] == "simulation":
            stats["simulation"] += 1
            continue
        if not event["reviewed"]:
            stats["unreviewed"] += 1
            continue
        r = copy.deepcopy(event["sample"])
        key = (r["seq"], event["position_1based"])
        if key in used:
            raise ValueError(f"Repeated sequence/site needs manual reconciliation, not extra weight: {key}")
        used.add(key)
        mask = [int(i + 1 == event["position_1based"]) for i in range(len(r["seq"]))]
        r.update(record_id="feedback_" + event["observation_id"], label_mask=mask,
                 tested_mask=mask.copy(), ca_label=[event["outcome"] if m else 0 for m in mask],
                 confirmed_negative_mask=[int(m and event["outcome"] == 0) for m in mask],
                 experiment_confirmed=True, evidence_type="wet_experiment", supervision_weight=1.0,
                 parent_group=event["parent_group"], parent_start_1based=None,
                 feedback_observation=event)
        rows.append(r)
        stats["accepted"] += 1
    if held_groups & {r["parent_group"] for r in rows}:
        raise ValueError("Holdout parent overlaps training; curate group split before retraining")
    return rows, dict(stats)


def evaluate(y, scores, passes=None):
    result = core.metrics(y, scores, passes)
    if len(set(y)) < 2:
        for key in ("accuracy", "precision", "f1", "average_precision"):
            result[key] = None
    return result


def train(args):
    if args.epochs < 1 or args.output.exists():
        raise ValueError("Positive epochs and a new output directory required")
    rows, feedback_stats = merge_feedback(prepare(args.source, args.parents), args.feedback)
    torch.set_num_threads(4)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    test_sequence = "AGCKNFFWKTFTSC"
    test_overlap = [r["record_id"] for r in rows
                    if r["seq"] in test_sequence or test_sequence in r["seq"]]
    config = {"version": "zqy8-v1",
              "source_directory": relative_record(args.source),
              "test_2mi1_sequence_overlap_records": test_overlap,
              "seed": args.seed, "epochs": args.epochs,
              "hard_rules": RULES, "model": "frozen official ProteinMPNN + 363-input linear site head",
              "encoder_device": str(device), "torch_version": str(torch.__version__),
              "gpu": torch.cuda.get_device_name(0) if device.type == "cuda" else None,
              "source_sha256": {name: core.sha(args.source / name) for name in (positive_filename(args.source), "训练集负例.jsonl")},
              "pretrained_sha256": core.sha(args.pretrained), "parents_sha256": core.sha(args.parents),
              "code_sha256": {p.name: core.sha(p) for p in (Path(__file__), Path(core.__file__))},
              "feedback_sha256": core.sha(args.feedback) if args.feedback else None,
              "weights": {"main_positive": 1.0, "confirmed_negative": 1.0, "low_precision": 0,
                          "validation": 0, "simulation": 0},
              "optimizer": {"name": "AdamW", "lr": 0.01, "weight_decay": 0.01},
              "selection": "fixed final epoch and 0.5 threshold; no validation tuning",
              "iteration": "full replay from frozen official encoder; never train on simulations/holdout",
              "efficiency": "archived metadata only; binary head is not a yield regressor",
              "label_semantics": "base source ca_label is tested mask; file provenance determines outcome",
              "limitations": ["Not an insertion generator", "Main positive file authoritative; full-atom feature provenance not independently audited",
                              "Negative parents concentrated; homology not audited", "2MI1 is not evaluated: reaction-state structure unresolved; known ordering must not guide tuning",
                              "Pro exclusion is user-adopted hypothesis", "nearby_atoms is learned input, not an invented target-range penalty"]}
    args.output.mkdir(parents=True)
    core.write_json(args.output / "config.json", config)
    core.write_jsonl(args.output / "prepared/training.jsonl", rows)
    if args.feedback:
        core.write_jsonl(args.output / "prepared/feedback_archive.jsonl", core.load_jsonl(args.feedback))
    backbone = core.make_backbone(args.pretrained, device)
    vectors, targets, weights, sites = [], [], [], []
    for r in rows:
        f = core.features(backbone, r, device)
        for i, mask in enumerate(r["label_mask"]):
            if not mask:
                continue
            vectors.append(f[i]); targets.append(r["ca_label"][i]); weights.append(r["supervision_weight"])
            sites.append({"record_id": r["record_id"], "seq": r["seq"], "position_1based": i+1,
                          "parent_group": r["parent_group"], "y": r["ca_label"][i],
                          "rejection_reasons": reasons(r, i, RULES), **{k: r[k][i] for k in core.BIO}})
    x = torch.stack(vectors); y = torch.tensor(targets, dtype=torch.float32); w = torch.tensor(weights)
    if len(set(targets)) != 2:
        raise ValueError("Two experimental outcome classes required")
    columns = list(range(x.shape[1]))
    head, logs = core.fit(x, y, w, args.epochs, columns, args.seed)
    scores = core.score_head(head, x, columns)
    passes = np.array([not s["rejection_reasons"] for s in sites])
    training = {"raw": evaluate(targets, scores), "filtered": evaluate(targets, scores, passes),
                "rules_only": evaluate(targets, passes.astype(float)), "scope": "training resubstitution only"}
    groups = np.array([s["parent_group"] for s in sites])
    diagnostics = {}
    for name, cols in {"full": columns, "biochemical_only": list(range(359, 363)),
                       "without_biochemical": list(range(359))}.items():
        folds, held_predictions = [], []
        for group in sorted(set(groups)):
            tr = np.where(groups != group)[0]; te = np.where(groups == group)[0]
            fold = {"parent": str(group), "train_n": len(tr), "held_n": len(te)}
            if len(set(targets[i] for i in tr)) < 2:
                folds.append({**fold, "status": "skipped_single_class_training"})
                continue
            model, _ = core.fit(x[tr], y[tr], w[tr], args.epochs, cols, args.seed)
            pred = core.score_head(model, x[te], cols)
            folds.append({**fold, "status": "evaluated", "raw": evaluate(y[te].tolist(), pred),
                          "filtered": evaluate(y[te].tolist(), pred, passes[te])})
            held_predictions.extend({**sites[i], "raw_score": float(p)} for i, p in zip(te, pred))
        diagnostics[name] = {"folds": folds, "pooled": evaluate([s["y"] for s in held_predictions], [s["raw_score"] for s in held_predictions]),
                             "scope": "grouped diagnostic; single-class test cannot estimate false positives"}
        core.write_jsonl(args.output / f"diagnostics/{name}_heldout.jsonl", held_predictions)
    checkpoint = args.output / "final_model/site_predictor.pt"
    checkpoint.parent.mkdir()
    torch.save({"backbone_state_dict": {k: v.cpu() for k, v in backbone.state_dict().items()},
                "head_state_dict": head.state_dict(), "columns": columns, "config": config}, checkpoint)
    saved = torch.load(checkpoint, weights_only=True, map_location="cpu")
    reload_head = torch.nn.Linear(len(columns), 1)
    reload_head.load_state_dict(saved["head_state_dict"])
    if not np.array_equal(scores, core.score_head(reload_head, x, columns)):
        raise AssertionError("Reload mismatch")
    core.write_jsonl(args.output / "training_predictions.jsonl", [{**s, "raw_score": float(p), "passes_hard_filter": bool(ok)} for s, p, ok in zip(sites, scores, passes)])
    core.write_json(args.output / "training_log.json", logs)
    core.write_json(args.output / "diagnostics/grouped_ablation.json", diagnostics)
    summary = {"training_records": len(rows), "training_sites": len(sites), "positive_sites": sum(targets),
               "negative_sites": len(targets)-sum(targets), "parent_groups": dict(Counter(groups)),
               "trainable_parameters": sum(p.numel() for p in head.parameters()),
               "training": training, "diagnostics": diagnostics, "feedback": feedback_stats,
               "checkpoint_sha256": core.sha(checkpoint), "reload_verified": True}
    core.write_json(args.output / "summary.json", summary)
    for name, digest in config["source_sha256"].items():
        if core.sha(args.source / name) != digest:
            raise AssertionError("Source changed during training")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


def predict(args):
    if args.output.exists():
        raise FileExistsError(args.output)
    torch.set_num_threads(4)
    saved = torch.load(args.checkpoint, weights_only=True, map_location="cpu")
    if saved["config"].get("version") not in ("zqy7-v1", "zqy8-v1"):
        raise ValueError("Use a ZQY7/8 checkpoint; do not silently change old model rules")
    backbone = core.ProteinMPNN(num_letters=21, node_features=128, edge_features=128, hidden_dim=128,
                               num_encoder_layers=3, num_decoder_layers=3, vocab=21, k_neighbors=48,
                               augment_eps=0.0, dropout=0.0, ca_only=True).eval().requires_grad_(False)
    backbone.load_state_dict(saved["backbone_state_dict"], strict=True)
    head = torch.nn.Linear(len(saved["columns"]), 1)
    head.load_state_dict(saved["head_state_dict"])
    output = []
    for row in core.load_jsonl(args.input):
        core.validate(row)
        mask = row.get("feature_mask", row.get("tested_mask"))
        if mask is None or len(mask) != len(row["seq"]) or any(v not in (0, 1) for v in mask):
            raise ValueError("Explicit feature_mask required; zero SASA may be a valid measurement")
        score = core.score_head(head, core.features(backbone, row, torch.device("cpu")), saved["columns"])
        for i, aa in enumerate(row["seq"]):
            if aa == "C" and mask[i]:
                why = reasons(row, i, saved["config"]["hard_rules"])
                output.append({"record_id": row.get("record_id"), "seq": row["seq"], "position_1based": i+1,
                               "raw_score": float(score[i]), "passes_hard_filter": not why,
                               "predicted_positive": bool(score[i] >= 0.5 and not why), "rejection_reasons": why})
    core.write_jsonl(args.output, output)
    print(f"Scored {len(output)} existing Cys sites; scores are not calibrated efficiencies.")



def relative_record(path):
    """Store portable provenance for files within the repository."""
    import os
    return Path(os.path.relpath(Path(path).resolve(), ROOT)).as_posix()
