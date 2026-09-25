"""Cysteine insertion design with explicit candidate-to-parent position mapping."""
import hashlib
import json
import re
import csv
import tempfile
import math
from pathlib import Path
from types import SimpleNamespace

AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")


def validate_sequence(sequence):
    if not isinstance(sequence, str) or len(sequence) < 2:
        raise ValueError("Input peptide must contain at least two standard amino acids")
    if any(aa not in AMINO_ACIDS for aa in sequence):
        raise ValueError("Use uppercase one-letter codes for the 20 standard amino acids")
    return sequence


def insertion_candidates(sequence, record_id="peptide"):
    """Enumerate real insertions: parent length N becomes N+1, never substitution.

    Keep separate site mappings for insertions beside existing Cys, even when
    candidate sequences coincide. Ranking must deduplicate final sequences.
    """
    validate_sequence(sequence)
    if not isinstance(record_id, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,79}", record_id):
        raise ValueError("record_id must be 1-80 safe letters/digits/underscore/dot/hyphen")
    rows = []
    for gap in range(len(sequence) + 1):
        candidate = sequence[:gap] + "C" + sequence[gap:]
        rows.append({
            "candidate_id": f"{record_id}_insert_C_after_{gap}",
            "parent_id": record_id, "parent_sequence": sequence,
            "sequence": candidate, "insert_after_parent_position_1based": gap,
            "inserted_position_1based": gap + 1,
            "candidate_to_parent_1based": list(range(1, gap + 1)) + [None] + list(range(gap + 1, len(sequence) + 1)),
            "sequence_sha256": hashlib.sha256(candidate.encode("ascii")).hexdigest(),
            "previous_residue": sequence[gap - 1] if gap else None,
            "n_terminal_insertion": gap == 0,
            "status": "awaiting_candidate_structure",
        })
    return rows


def rank_candidates(rows, top_k=3):
    """Prefer passing candidates, then fill with scored warnings; never invent scores."""
    if top_k < 1:
        raise ValueError("top-k must be positive")
    eligible = sorted((r for r in rows if isinstance(r.get("prediction_score"), (int, float))
                       and math.isfinite(r["prediction_score"])),
                      key=lambda r: (not r["eligible"], -r["prediction_score"], r["candidate_id"]))
    selected, seen = [], set()
    for row in eligible:
        if row["sequence"] in seen:
            continue
        seen.add(row["sequence"])
        selected.append({**row, "rank": len(selected) + 1,
                         "selection_status": "recommended" if row["eligible"] else "fallback",
                         "warning": candidate_warning(row)})
        if len(selected) == top_k:
            break
    return selected


def candidate_warning(row):
    if row["eligible"]:
        return ""
    meanings = {
        "free_N_terminal_Cys_has_no_preceding_backbone_amide": "新增 Cys 位于游离 N 端，缺少前接主链酰胺",
        "structure_geometry_check_failed": "预测结构未通过几何检查，需修复并重新评估",
        "disulfide": "Cys 可能参与二硫键",
        "previous_residue_proline": "前位为 Pro，触发现有保守过滤规则",
    }
    reasons = "；".join(meanings.get(r, r) for r in row.get("exclusion_reasons", []))
    return "补足候选：根据已有了解，该序列的成功可能性较低，尚未经实验验证。原因：" + (reasons or "未通过既有过滤规则")


def score(args):
    from . import core, pipeline, barrier, structure
    if args.output.exists():
        raise FileExistsError(args.output)
    if args.top_k < 1:
        raise ValueError("top-k must be positive")
    candidates = core.load_jsonl(args.candidates)
    if not candidates:
        raise ValueError("Empty candidate list")
    # Recompute mapping, rather than trusting an edited position/sequence pairing.
    expected = insertion_candidates(candidates[0]["parent_sequence"], candidates[0]["parent_id"])
    if candidates != expected:
        raise ValueError("Candidate manifest changed; regenerate with design.py prepare")
    manifest = json.loads(args.structures.read_text(encoding="utf-8-sig"))
    if not isinstance(manifest, list) or not manifest:
        raise ValueError("Structure manifest must be a nonempty list")
    by_seq = {}
    for entry in manifest:
        seq = entry.get("sequence")
        if seq not in {r["sequence"] for r in candidates} or seq in by_seq:
            raise ValueError("Unknown or duplicate structure-manifest sequence")
        for key in ("structure_file", "chain_id", "structure_method", "structure_notes"):
            if not isinstance(entry.get(key), str) or not entry[key].strip():
                raise ValueError(f"Structure entry requires {key}")
        path = Path(entry["structure_file"])
        # Structure references are relative to their manifest, not current cwd.
        path = path if path.is_absolute() else args.structures.parent / path
        if not path.is_file():
            raise FileNotFoundError(path)
        by_seq[seq] = {**entry, "path": path.resolve()}
    missing = [r["candidate_id"] for r in candidates if r["sequence"] not in by_seq]
    if missing and not args.allow_partial:
        raise ValueError("Missing candidate structures; supply all or explicitly use --allow-partial")
    samples, evidence_rows, chains = [], {}, {}
    for candidate in candidates:
        seq = candidate["sequence"]
        if seq not in by_seq:
            continue
        entry = by_seq[seq]
        sample, evidence, chain = structure.features(entry["path"], seq, candidate["inserted_position_1based"],
                                                     entry["chain_id"], entry.get("model_index", 1))
        sample.update(record_id=candidate["candidate_id"], parent_group=candidate["parent_id"])
        samples.append(sample)
        evidence.update(source_structure_sha256=core.sha(entry["path"]),
                        structure_method=entry["structure_method"], structure_notes=entry["structure_notes"])
        evidence_rows[candidate["candidate_id"]] = evidence
        chains[candidate["candidate_id"]] = chain
    if not samples:
        raise ValueError("No candidate structures available")
    with tempfile.TemporaryDirectory() as temp:
        inputs, outputs = Path(temp) / "input.jsonl", Path(temp) / "scores.jsonl"
        core.write_jsonl(inputs, samples)
        pipeline.predict(SimpleNamespace(input=inputs, checkpoint=args.checkpoint, output=outputs))
        scored = {r["record_id"]: r for r in core.load_jsonl(outputs)}
    full = []
    sample_by_id = {r["record_id"]: r for r in samples}
    for candidate in candidates:
        identity = candidate["candidate_id"]
        if identity not in scored:
            full.append({**candidate, "eligible": False, "prediction_score": None,
                         "exclusion_reasons": ["missing_candidate_structure"]})
            continue
        s, e = scored[identity], evidence_rows[identity]
        reasons = list(s["rejection_reasons"])
        if candidate["n_terminal_insertion"]:
            reasons.append("free_N_terminal_Cys_has_no_preceding_backbone_amide")
        if not e["geometry_passed"]:
            reasons.append("structure_geometry_check_failed")
        full.append({**candidate, "status": "scored", "prediction_score": s["raw_score"],
                     "classified_positive_at_0_5": s["predicted_positive"],
                     "eligible": not reasons, "exclusion_reasons": reasons, "structure_evidence": e})
    args.output.mkdir(parents=True)
    (args.output / "structures").mkdir()
    for row in full:
        identity = row["candidate_id"]
        if identity not in chains:
            continue
        destination = args.output / "structures" / f"{identity}.pdb"
        structure.write_chain(chains[identity], destination)
        row["structure_file"] = f"structures/{identity}.pdb"
        sample_by_id[identity]["structure_file"] = str(destination.resolve())
    # Re-select to retain the copied structure paths added above.
    top = rank_candidates(full, args.top_k)
    fields = ["rank", "candidate_id", "sequence", "insert_after_parent_position_1based",
              "inserted_position_1based", "prediction_score", "classified_positive_at_0_5", "structure_file",
              "eligible", "selection_status", "warning", "exclusion_reasons"]
    with (args.output / "top3.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(top)
    core.write_jsonl(args.output / "top3.jsonl", top)
    core.write_jsonl(args.output / "all_candidates.jsonl", full)
    core.write_jsonl(args.output / "model_inputs.jsonl", samples)
    top_samples = [sample_by_id[r["candidate_id"]] for r in top]
    core.write_jsonl(args.output / "top3.model_inputs.jsonl", top_samples)
    if top_samples:
        barrier.requests(SimpleNamespace(input=args.output / "top3.model_inputs.jsonl", output=args.output / "barrier_requests.jsonl"))
    core.write_json(args.output / "run.json", {
        "operation": "single_Cys_insertion_design", "candidate_count": len(candidates),
        "scored_count": len(scored), "unique_top_sequences": len(top), "requested_top_k": args.top_k,
        "missing_candidate_ids": missing, "ranking_scope": "available structures only" if missing else "all enumerated insertion sites",
        "model_sha256": core.sha(args.checkpoint), "candidate_manifest_sha256": core.sha(args.candidates),
        "code_sha256": {"design": core.sha(__file__), "structure": core.sha(structure.__file__)},
        "structure_manifest_sha256": core.sha(args.structures),
        "fallback_count": sum(r["selection_status"] == "fallback" for r in top),
        "ranking": "passing candidates first, then scored filtered candidates to fill top-k; descending score within each tier; unique sequences",
        "threshold_policy": "0.5 classification reported separately; not a threshold for design prioritization",
        "limitations": ["No demonstrated insertion-design generalization", "Scores are not cyclization efficiency",
                        "Structure-feature calculation may differ from historical training data",
                        "Geometry checks do not validate folding or reaction chemistry",
                        "Legacy previous-Pro hard filter retained from checkpoint; hypothesis only"]})
    report = ["# Cys 插入候选结构报告", "", f"枚举 {len(candidates)} 个插入位置，评分 {len(scored)} 个，输出 {len(top)} 条不同序列。", "",
              "优先选择通过过滤的候选；不足时按内部评分补足并提示风险。各组内按分数降序，序列去重。模型分数不是成功概率。", "",
              "结构检查仅含重原子完整性、主链键长与严重非局部重叠，不能证明折叠正确或反应可行。历史训练特征口径不完整，存在分布偏移风险。", ""]
    for r in top:
        e = r["structure_evidence"]
        report += [f"## 候选 {r['rank']}", "", f"序列：`{r['sequence']}`",
                   f"在原序列第 {r['insert_after_parent_position_1based']} 位之后插入 C；新位点编号 {r['inserted_position_1based']}。",
                   f"模型分数：{r['prediction_score']:.6f}；0.5 阈值分类：{r['classified_positive_at_0_5']}。",
                   f"Cys 侧链 SASA：{e['sasa_sidechain_A2']:.3f} Å²；SG 周围 5 Å 非氢原子数：{e['nearby_atoms_5A']}。",
                   f"结构来源：{e['structure_method']}。结构文件：{r['structure_file']}。", ""]
        if r["warning"]:
            report += [f"**{r['warning']}**", ""]
    if len(top) < args.top_k:
        report += ["已有分数的不同序列不足请求数量；已用过滤候选补足至可用数量，不重复序列或虚构缺失分数。", ""]
    if missing:
        report += ["部分候选没有结构，本次排名只覆盖已提供结构的候选。", ""]
    (args.output / "report.md").write_text("\n".join(report), encoding="utf-8")
    print(f"Exported {len(top)} distinct insertion designs to {args.output}")


def prepare(sequence, record_id, output):
    if output.exists():
        raise FileExistsError(output)
    rows = insertion_candidates(sequence, record_id)
    output.mkdir(parents=True)
    (output / "candidates.jsonl").write_text(
        "".join(json.dumps(r, ensure_ascii=False) + "\n" for r in rows), encoding="utf-8")
    # Identical sequences need only one structure prediction, but all insertion
    # interpretations remain represented in candidates.jsonl.
    unique = {}
    for r in rows:
        unique.setdefault(r["sequence"], r["candidate_id"])
    (output / "candidate_sequences.fasta").write_text(
        "".join(f">{identity}\n{seq}\n" for seq, identity in unique.items()), encoding="utf-8")
    manifest = [{"sequence": seq, "structure_id": identity,
                 "candidate_ids": [r["candidate_id"] for r in rows if r["sequence"] == seq],
                 "structure_file": None, "chain_id": None, "model_index": 1,
                 "structure_method": None, "structure_notes": None}
                for seq, identity in unique.items()]
    (output / "structure_manifest.template.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    (output / "preparation.json").write_text(json.dumps({
        "parent_id": record_id, "parent_sequence": sequence,
        "operation": "insert_one_C", "candidate_count": len(rows),
        "unique_sequence_count": len(unique),
        "status": "awaiting_candidate_structures; no model scores or Top3 assigned",
    }, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return rows
