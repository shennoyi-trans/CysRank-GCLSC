"""External simulation contract and a separate, protocol-specific barrier regressor.

No quantum chemistry is run here. Simulated barriers never become experimental
binary labels. The existing frozen encoder supplies label-blind site features.
"""
import copy
import hashlib
import json
import math
from collections import Counter

import numpy as np
import torch

from . import core

VERSION = "cysrank-barrier-v1"
FACTORS = {"kcal/mol": 1.0, "kJ/mol": 1 / 4.184,
           "eV": 23.0605478306, "Hartree": 627.5094740631}
PROTOCOL_FIELDS = ("software", "software_version", "method", "reaction_id",
                   "reaction_step", "reference_state", "standard_state", "solvent",
                   "protonation_policy", "conformer_policy", "settings_reference")


def nonempty(value, name):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"Nonempty {name} required")


def finite(value, name):
    if type(value) not in (float, int) or not math.isfinite(value):
        raise ValueError(f"Finite numeric {name} required")


def digest(value):
    return hashlib.sha256(json.dumps(value, ensure_ascii=False, sort_keys=True,
                                    separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def sample_digest(sample):
    return digest({k: sample[k] for k in ("seq", "xyz", *core.BIO)})


def validate_barrier(event):
    """Called after the common feedback validator; returns a normalized copy."""
    if event["evidence_type"] != "simulation" or event.get("outcome") is not None:
        raise ValueError("Barrier observations require simulation evidence and outcome=null")
    b = event["barrier"]
    if not isinstance(b, dict) or b.get("schema_version") != VERSION:
        raise ValueError(f"barrier.schema_version must be {VERSION}")
    if b.get("status") not in ("completed", "failed"):
        raise ValueError("barrier.status must be completed or failed; pending requests cannot be ingested")
    if b.get("quantity") not in ("delta_e_dagger", "delta_g_dagger"):
        raise ValueError("Specify delta_e_dagger or delta_g_dagger, not a total/reaction energy")
    if b.get("unit") not in FACTORS:
        raise ValueError(f"Supported barrier units: {list(FACTORS)}")
    if type(b.get("quality_passed")) is not bool:
        raise ValueError("Explicit barrier.quality_passed boolean required")
    nonempty(b.get("quality_notes"), "barrier.quality_notes")
    nonempty(b.get("structure_reference"), "barrier.structure_reference")
    if b.get("sample_sha256") != sample_digest(event["sample"]):
        raise ValueError("sample_sha256 does not match the supplied model-input structure/features")
    protocol = b.get("protocol")
    if not isinstance(protocol, dict):
        raise ValueError("barrier.protocol object required")
    for key in PROTOCOL_FIELDS:
        nonempty(protocol.get(key), f"protocol.{key}")
    finite(protocol.get("temperature_kelvin"), "protocol.temperature_kelvin")
    if protocol["temperature_kelvin"] <= 0:
        raise ValueError("Temperature must be positive Kelvin")
    if type(b.get("charge")) is not int or type(b.get("multiplicity")) is not int or b["multiplicity"] < 1:
        raise ValueError("Integer charge and positive integer multiplicity required")
    if b["status"] == "failed":
        if b.get("value") is not None or b["quality_passed"]:
            raise ValueError("Failed calculations require value=null and quality_passed=false")
        nonempty(b.get("failure_reason"), "barrier.failure_reason")
        return None
    finite(b.get("value"), "barrier.value")
    value = b["value"] * FACTORS[b["unit"]]
    finite(value, "converted barrier")
    # Signed values are retained: do not silently clamp or relabel them.
    return value


def protocol_key(event):
    b = event["barrier"]
    return digest({"quantity": b["quantity"], "protocol": b["protocol"]})


def eligible(event):
    b = event.get("barrier", {})
    return (event.get("evidence_type") == "simulation" and event.get("reviewed") is True
            and b.get("status") == "completed" and b.get("quality_passed") is True)


def collect(events):
    from .pipeline import validate_event
    seen, rows = set(), []
    stats = Counter(total=len(events), training=0, holdout=0, ineligible=0)
    # Validate every input, even excluded records. Check split leakage before filtering.
    for e in events:
        validate_event(e)
        if e["observation_id"] in seen:
            raise ValueError("Duplicate observation_id")
        seen.add(e["observation_id"])
        if "barrier" not in e:
            raise ValueError("Barrier training requires barrier observations only")
    groups = {e["parent_group"] for e in events if e["split"] == "train"}
    held = [e for e in events if e["split"] == "holdout"]
    if groups & {e["parent_group"] for e in held}:
        raise ValueError("Holdout parent overlaps training")
    train_seq = {e["sample"]["seq"] for e in events if e["split"] == "train"}
    if any(a in e["sample"]["seq"] or e["sample"]["seq"] in a for e in held for a in train_seq):
        raise ValueError("Holdout sequence overlaps training; curate the split")
    used = set()
    for e in events:
        if not eligible(e):
            stats["ineligible"] += 1
            continue
        # One reviewed aggregate per sequence/site/protocol prevents repeated runs
        # or conformer trials from silently increasing a site's weight.
        key = (e["sample"]["seq"], e["position_1based"], protocol_key(e))
        if key in used:
            raise ValueError("Repeated sequence/site/protocol; aggregate or reconcile simulations first")
        used.add(key)
        rows.append(e)
        stats["training" if e["split"] == "train" else "holdout"] += 1
    if len({protocol_key(e) for e in rows}) > 1:
        raise ValueError("Mixed quantities or protocols: train separate barrier models")
    return rows, dict(stats)


def regression_metrics(y, prediction):
    if not len(y):
        return None
    error = np.asarray(prediction) - np.asarray(y)
    return {"n": len(y), "mae_kcal_mol": float(np.abs(error).mean()),
            "rmse_kcal_mol": float(np.sqrt(np.square(error).mean()))}


def fit_ridge(x, y, alpha):
    finite(alpha, "alpha")
    if alpha <= 0:
        raise ValueError("alpha must be positive")
    # Fit scaling on train only. Penalize weights, not the intercept.
    mean = x.mean(axis=0)
    scale = x.std(axis=0)
    scale[scale < 1e-8] = 1.0
    z = (x - mean) / scale
    center = float(y.mean())
    weights = np.linalg.solve(z.T @ z + alpha * np.eye(x.shape[1]), z.T @ (y - center))
    return {"mean": mean.tolist(), "scale": scale.tolist(),
            "weights": weights.tolist(), "intercept": center}


def apply_regressor(model, x):
    return ((x - np.asarray(model["mean"])) / np.asarray(model["scale"])) @ np.asarray(model["weights"]) + model["intercept"]


def vectors(rows, pretrained):
    torch.set_num_threads(4)
    device = torch.device("cpu")
    backbone = core.make_backbone(pretrained, device)
    return np.stack([core.features(backbone, e["sample"], device)[e["position_1based"] - 1].numpy()
                     for e in rows]).astype(np.float64)


def train(args):
    if args.output.exists():
        raise FileExistsError(args.output)
    rows, stats = collect(core.load_jsonl(args.input))
    training = [e for e in rows if e["split"] == "train"]
    if len(training) < 2:
        raise ValueError("At least two distinct reviewed, quality-passed training sites required")
    targets = np.array([validate_barrier(e) for e in rows])
    mask = np.array([e["split"] == "train" for e in rows])
    x = vectors(rows, args.pretrained)
    head = fit_ridge(x[mask], targets[mask], args.alpha)
    model = {"version": VERSION, "kind": "frozen_ProteinMPNN_ridge_barrier",
             "quantity": rows[0]["barrier"]["quantity"], "unit": "kcal/mol",
             "protocol": rows[0]["barrier"]["protocol"], "protocol_sha256": protocol_key(rows[0]),
             "pretrained_sha256": core.sha(args.pretrained), "alpha": args.alpha,
             "head": head, "feedback_sha256": core.sha(args.input),
             "code_sha256": {"barrier": core.sha(__file__), "features": core.sha(core.__file__)},
             "training_observations": [e["observation_id"] for e in training],
             "training_parent_groups": sorted({e["parent_group"] for e in training}),
             "scope": "simulation surrogate only; not experimental probability or efficiency"}
    prediction = apply_regressor(head, x)
    summary = {"feedback": stats, "training": regression_metrics(targets[mask], prediction[mask]),
               "holdout": regression_metrics(targets[~mask], prediction[~mask]),
               "unit": "kcal/mol", "scope": model["scope"],
               "validation_note": "No independent holdout supplied" if mask.all() else
               "Holdout excluded from fitting and normalization; homology needs external review"}
    args.output.mkdir(parents=True)
    core.write_json(args.output / "barrier_model.json", model)
    core.write_json(args.output / "protocol.json", {"quantity": model["quantity"], "protocol": model["protocol"]})
    core.write_json(args.output / "summary.json", summary)
    core.write_jsonl(args.output / "accepted_feedback.jsonl", rows)
    core.write_jsonl(args.output / "predictions.jsonl", [
        {"observation_id": e["observation_id"], "split": e["split"],
         "target_kcal_mol": float(t), "predicted_barrier_kcal_mol": float(p)}
        for e, t, p in zip(rows, targets, prediction)])
    print(json.dumps(summary, ensure_ascii=False, indent=2))


def predict(args):
    if args.output.exists():
        raise FileExistsError(args.output)
    model = json.loads(args.model.read_text(encoding="utf-8"))
    if model.get("version") != VERSION or model.get("kind") != "frozen_ProteinMPNN_ridge_barrier":
        raise ValueError("Unsupported barrier model")
    if core.sha(args.pretrained) != model["pretrained_sha256"]:
        raise ValueError("Pretrained encoder differs from the trained barrier model")
    # Explicit protocol declaration prevents accidentally ranking candidates under
    # a different method, reaction step, solvent or reference state.
    contract = json.loads(args.protocol.read_text(encoding="utf-8"))
    if digest(contract) != model["protocol_sha256"]:
        raise ValueError("Prediction quantity/protocol differs from the barrier model")
    rows = []
    for sample in core.load_jsonl(args.input):
        core.validate(sample)
        nonempty(sample.get("record_id"), "record_id")
        mask = sample.get("feature_mask")
        if not isinstance(mask, list) or len(mask) != len(sample["seq"]) or any(type(v) is not int or v not in (0, 1) for v in mask):
            raise ValueError("Explicit binary feature_mask required")
        for i, flag in enumerate(mask):
            if flag:
                if sample["seq"][i] != "C":
                    raise ValueError("Barrier prediction targets must be Cys")
                rows.append({"sample": sample, "position_1based": i + 1})
    if not rows:
        raise ValueError("No Cys sites selected")
    keys = [(e["sample"]["record_id"], e["position_1based"]) for e in rows]
    if len(set(keys)) != len(keys):
        raise ValueError("Duplicate prediction record/site")
    values = apply_regressor(model["head"], vectors(rows, args.pretrained))
    output = [{"record_id": e["sample"]["record_id"], "position_1based": e["position_1based"],
               "predicted_barrier_kcal_mol": float(p), "quantity": model["quantity"],
               "protocol_sha256": model["protocol_sha256"], "model_sha256": core.sha(args.model),
               "scope": model["scope"]} for e, p in zip(rows, values)]
    output.sort(key=lambda e: (e["predicted_barrier_kcal_mol"], e["record_id"], e["position_1based"]))
    core.write_jsonl(args.output, output)
    print(f"Predicted {len(output)} simulation barriers; ascending order within the supplied protocol.")


def requests(args):
    if args.output.exists():
        raise FileExistsError(args.output)
    output = []
    seen = set()
    for sample in core.load_jsonl(args.input):
        core.validate(sample)
        nonempty(sample.get("record_id"), "record_id")
        if sample["record_id"] in seen:
            raise ValueError("Duplicate record_id")
        seen.add(sample["record_id"])
        mask = sample.get("feature_mask")
        if not isinstance(mask, list) or len(mask) != len(sample["seq"]) or any(type(v) is not int or v not in (0, 1) for v in mask):
            raise ValueError("Explicit binary feature_mask required")
        snapshot = {k: copy.deepcopy(sample[k]) for k in ("seq", "xyz", *core.BIO)}
        for i, flag in enumerate(mask):
            if not flag:
                continue
            if sample["seq"][i] != "C":
                raise ValueError("Simulation request targets must be Cys")
            output.append({"observation_id": f"{sample['record_id']}:C{i+1}:{sample_digest(snapshot)[:12]}",
                           "evidence_type": "simulation", "split": None,
                           "parent_group": sample.get("parent_group"), "assay_id": None,
                           "conditions": {}, "provenance": None, "sample": snapshot,
                           "position_1based": i + 1, "outcome": None, "reviewed": False,
                           "barrier": {"schema_version": VERSION, "status": "pending",
                                       "value": None, "unit": "kcal/mol", "quantity": None,
                                       "quality_passed": False, "quality_notes": None,
                                       "sample_sha256": sample_digest(snapshot),
                                       "structure_reference": sample.get("structure_file"),
                                       "charge": None, "multiplicity": None,
                                       "protocol": {**{k: None for k in PROTOCOL_FIELDS},
                                                    "temperature_kelvin": None}}})
    if not output:
        raise ValueError("No Cys sites selected")
    core.write_jsonl(args.output, output)
    print(f"Exported {len(output)} pending requests; fill metadata and results before ingestion.")
