"""Shared label-blind features, metrics and optimization for the ZQY8 site head."""
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import torch
import torch.nn.functional as F

from .protein_mpnn_utils import ProteinMPNN, gather_nodes

ROOT = Path(__file__).resolve().parents[1]
ALPHABET = "ACDEFGHIKLMNPQRSTVWYX"
BIO = ("is_disulfide", "sasa", "prev_is_pro", "nearby_atoms")


def load_jsonl(path):
    return [json.loads(line) for line in Path(path).read_text(encoding="utf-8-sig").splitlines() if line.strip()]


def write_json(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")


def write_jsonl(path, rows):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text("".join(json.dumps(r, ensure_ascii=False, allow_nan=False) + "\n" for r in rows), encoding="utf-8")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def validate(r):
    n = len(r["seq"])
    if n < 3 or any(c not in ALPHABET for c in r["seq"]):
        raise ValueError("Need at least three valid sequence residues for CA-only encoding")
    for key in ("xyz", *BIO):
        x = np.asarray(r[key], dtype=np.float32)
        shape = (n, 3) if key == "xyz" else (n,)
        if x.shape != shape or not np.isfinite(x).all():
            raise ValueError(f"Invalid {key}: expected {shape}, got {x.shape}")
        if key in ("is_disulfide", "prev_is_pro") and not np.isin(x, [0, 1]).all():
            raise ValueError(f"Invalid binary {key}")
        if key in ("sasa", "nearby_atoms") and (x < 0).any():
            raise ValueError(f"Negative {key}")
    if "ca_label" in r:
        x = np.asarray(r["ca_label"])
        if x.shape != (n,) or not np.isin(x, [0, 1]).all():
            raise ValueError("Invalid ca_label")


def make_backbone(weights, device):
    model = ProteinMPNN(num_letters=21, node_features=128, edge_features=128,
                        hidden_dim=128, num_encoder_layers=3, num_decoder_layers=3,
                        vocab=21, k_neighbors=48, augment_eps=0.0, dropout=0.0, ca_only=True)
    checkpoint = torch.load(weights, map_location="cpu", weights_only=True)
    model.load_state_dict(checkpoint["model_state_dict"], strict=True)
    model.to(device).eval().requires_grad_(False)
    return model


@torch.no_grad()
def features(model, r, device):
    """Never reads labels, provenance, sample weights or experimental outcomes."""
    validate(r)
    n = len(r["seq"])
    x = torch.tensor(r["xyz"], dtype=torch.float32, device=device)[None]
    mask = torch.ones(1, n, device=device)
    idx = torch.arange(n, device=device)[None]
    edges, neighbors = model.features(x, mask, idx, torch.zeros_like(idx))
    h = torch.zeros(1, n, 128, device=device)
    e = model.W_e(edges)
    attend = mask[:, :, None] * gather_nodes(mask[:, :, None], neighbors).squeeze(-1)
    for layer in model.encoder_layers:
        h, e = layer(h, e, neighbors, mask, attend)
    # Fixed per-residue normalization, no validation statistics fitted.
    h = F.layer_norm(h[0], (128,)).cpu()
    context = torch.zeros(n, 11, 21)
    for i in range(n):
        for offset in range(-5, 6):
            j = i + offset
            aa = r["seq"][j] if 0 <= j < n else "X"
            context[i, offset + 5, ALPHABET.index(aa)] = 1.0
    bio = torch.tensor([[r[k][i] for k in BIO] for i in range(n)], dtype=torch.float32)
    bio /= torch.tensor([1.0, 300.0, 1.0, 38.0])
    return torch.cat([h, context.flatten(1), bio], dim=1)


def metrics(y, scores, pass_mask=None):
    y = np.asarray(y, dtype=int)
    s = np.asarray(scores, dtype=float)
    pred = s >= 0.5
    if pass_mask is not None:
        pred &= np.asarray(pass_mask, dtype=bool)
    tp = int(np.sum(pred & (y == 1)))
    fp = int(np.sum(pred & (y == 0)))
    tn = int(np.sum(~pred & (y == 0)))
    fn = int(np.sum(~pred & (y == 1)))
    ap = None
    if len(set(y)) == 2 and pass_mask is None:
        # Average precision grouped by score threshold, with correct tie handling.
        recall_prev, ap = 0.0, 0.0
        for threshold in sorted(set(s), reverse=True):
            selected = s >= threshold
            true = int(np.sum(y[selected] == 1))
            recall = true / int(np.sum(y == 1))
            ap += (recall - recall_prev) * true / int(selected.sum())
            recall_prev = recall
    return {
        "n": len(y), "positive_sites": int(np.sum(y)), "negative_sites": int(np.sum(y == 0)),
        "tp": tp, "fp": fp, "tn": tn, "fn": fn,
        "accuracy": (tp + tn) / len(y) if len(y) else None,
        "precision": tp / (tp + fp) if tp + fp else None,
        "recall": tp / (tp + fn) if tp + fn else None,
        "specificity": tn / (tn + fp) if tn + fp else None,
        "f1": 2 * tp / (2 * tp + fp + fn) if 2 * tp + fp + fn else None,
        "average_precision": ap,
        "threshold": 0.5,
        "note": "Single-class sets cannot estimate both sensitivity and specificity; AP omitted."
                if len(set(y)) < 2 else "Average precision is threshold-grouped AP, not trapezoidal PR-AUC.",
    }


def joint_loss(head, x, y, weights, columns, ranking=None):
    logits = head(x[:, columns]).squeeze(-1)
    binary_sum = (F.binary_cross_entropy_with_logits(logits, y, reduction="none") * weights).sum()
    pair = logits.new_zeros(())
    if ranking is not None:
        # All conformers jointly represent one experimental observation.
        delta = head(ranking['preferred'][:, columns]) - head(ranking['other'][:, columns])
        pair = F.softplus(-delta).mean()
    observation_weight = 1.0 if ranking is not None else 0.0
    return (binary_sum + observation_weight * pair) / (weights.sum() + observation_weight), binary_sum / weights.sum(), pair


def fit(x, y, weights, epochs, columns, seed, ranking=None):
    torch.manual_seed(seed)
    head = torch.nn.Linear(len(columns), 1)
    torch.nn.init.zeros_(head.weight)
    torch.nn.init.zeros_(head.bias)
    opt = torch.optim.AdamW(head.parameters(), lr=0.01, weight_decay=0.01)
    logs = []
    for epoch in range(1, epochs + 1):
        opt.zero_grad()
        loss, _, _ = joint_loss(head, x, y, weights, columns, ranking)
        if not torch.isfinite(loss):
            raise RuntimeError("Non-finite training loss")
        loss.backward()
        opt.step()
        with torch.no_grad():
            score = head(x[:, columns]).squeeze(-1).sigmoid().numpy()
            total, post_loss, pair = joint_loss(head, x, y, weights, columns, ranking)
        logs.append({"epoch": epoch, "weighted_bce": float(post_loss), "joint_loss": float(total),
                     "pairwise_loss": float(pair), **metrics(y.numpy(), score)})
    return head.eval(), logs


def score_head(head, x, columns):
    with torch.no_grad():
        return head(x[:, columns]).squeeze(-1).sigmoid().numpy()
