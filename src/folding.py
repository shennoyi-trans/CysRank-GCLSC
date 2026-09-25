"""Local ESMFold adapter. No hosted inference or synthetic fallback coordinates."""
import json
from pathlib import Path

from .paths import ROOT

MODEL_ID = "facebook/esmfold_v1"
MODEL_REVISION = "75a3841ee059df2bf4d56688166c8fb459ddd97a"


class ESMFoldPredictor:
    def __init__(self, model_path=MODEL_ID, device="auto", chunk_size=32,
                 local_files_only=False, progress=print):
        import torch
        from transformers import EsmForProteinFolding

        if device not in ("auto", "cpu", "cuda"):
            raise ValueError("device must be auto, cpu or cuda")
        if chunk_size < 1:
            raise ValueError("chunk-size must be positive")
        # The 3B language-model stem alone exceeds the comfortable budget of a 6 GB card.
        if device == "auto":
            device = "cuda" if torch.cuda.is_available() and torch.cuda.mem_get_info()[0] >= 12 * 1024**3 else "cpu"
        if device == "cuda" and not torch.cuda.is_available():
            raise RuntimeError("CUDA unavailable; use --device cpu")
        progress(f"加载 ESMFold（{device}）；首次使用需要下载模型权重，后续复用本地缓存。")
        self.model = EsmForProteinFolding.from_pretrained(
            str(model_path), cache_dir=str(ROOT / "models" / "folding-cache"),
            revision=MODEL_REVISION if str(model_path) == MODEL_ID else None,
            local_files_only=local_files_only, low_cpu_mem_usage=True,
            # Official v1 weights are a .bin checkpoint. Avoid Hub's background
            # safetensors conversion request, especially for offline operation.
            use_safetensors=False if str(model_path) == MODEL_ID else None,
        ).eval()
        # Keep the geometry trunk in float32; use the pretrained stem's half precision.
        self.model.esm.half()
        self.model.trunk.set_chunk_size(chunk_size)
        self.model.to(device)
        import transformers
        self.metadata = {"model": str(model_path), "revision": getattr(self.model.config, "_commit_hash", None),
                         "device": device, "chunk_size": chunk_size,
                         "torch_version": torch.__version__, "transformers_version": transformers.__version__,
                         "stem_dtype": str(next(self.model.esm.parameters()).dtype),
                         "method": "ESMFold v1 / Hugging Face Transformers"}

    def predict(self, sequence):
        import torch
        with torch.inference_mode():
            return self.model.infer_pdb(sequence)


def predict_manifest(preparation, predictor, progress=print):
    """Predict each unique sequence once and validate residue/atom correspondence."""
    from .structure import read_peptide
    entries = json.loads((preparation / "structure_manifest.template.json").read_text(encoding="utf-8"))
    folder = preparation / "predicted_structures"
    folder.mkdir()
    for index, entry in enumerate(entries, 1):
        progress(f"结构预测 {index}/{len(entries)}：{entry['sequence']}")
        pdb = predictor.predict(entry["sequence"])
        if not isinstance(pdb, str) or "ATOM" not in pdb:
            raise ValueError("Structure model returned no PDB atoms")
        path = folder / f"{entry['structure_id']}.pdb"
        path.write_text(pdb, encoding="utf-8")
        read_peptide(path, entry["sequence"], "A")
        entry.update(structure_file=path.relative_to(preparation).as_posix(), chain_id="A",
                     structure_method=predictor.metadata["method"],
                     structure_notes="Local single-sequence prediction; unrelaxed single conformer; free termini; no explicit reaction reagent. "
                                     + json.dumps(predictor.metadata, ensure_ascii=False))
    manifest = preparation / "structures.json"
    manifest.write_text(json.dumps(entries, ensure_ascii=False, indent=2), encoding="utf-8")
    return manifest


def run_design(args, progress=print):
    """Sequence -> structure -> original checkpoint -> self-contained external handoff."""
    import gc
    import shutil
    from types import SimpleNamespace
    from . import design

    if args.output.exists():
        raise FileExistsError(args.output)
    design.validate_sequence(args.sequence)
    if args.top_k < 1:
        raise ValueError("top-k must be positive")
    progress("枚举所有 Cys 插入位置")
    preparation = args.output / "preparation"
    design.prepare(args.sequence, args.record_id, preparation)
    predictor = ESMFoldPredictor(args.fold_model, args.device, args.chunk_size,
                                 args.local_files_only, progress)
    manifest = predict_manifest(preparation, predictor, progress)
    metadata = predictor.metadata
    del predictor
    gc.collect()
    import torch
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
    progress("提取结构特征并使用 CysRank 评分")
    output = args.output / "results"
    design.score(SimpleNamespace(candidates=preparation / "candidates.jsonl", structures=manifest,
                                checkpoint=args.checkpoint, output=output, top_k=args.top_k,
                                allow_partial=False))
    # Export a portable copy; keep original requests intact for their input hashes.
    requests = output / "barrier_requests.jsonl"
    if requests.exists():
        rows = [json.loads(line) for line in requests.read_text(encoding="utf-8").splitlines() if line]
        for row in rows:
            # sample_sha256 remains bound to the original model input, not a rewritten sample.
            row["barrier"]["structure_reference"] = "structures/" + Path(row["barrier"]["structure_reference"]).name
        (output / "external_barrier_requests.jsonl").write_text(
            "".join(json.dumps(row, ensure_ascii=False) + "\n" for row in rows), encoding="utf-8")
    (output / "folding.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
    (output / "EXTERNAL_EVALUATION.md").write_text(
        "# 外部 GFN2-xTB 评估交接\n\n"
        "top3.csv / top3.jsonl 为模型排序；warning 标记补足候选及原因。分数不是实验成功概率。\n\n"
        "structures/ 为预测肽的重原子结构。external_barrier_requests.jsonl 的 structure_reference 相对本文件；"
        "sample 保留原始模型输入用于哈希追溯。不要修改 sample 或 sample_sha256。\n\n"
        "计算方需补氢、设定端基/质子化、电荷、溶剂、试剂及明确反应步骤，构建反应物/产物及路径，"
        "核验过渡态后报告 ΔE‡ 或 ΔG‡。本项目未执行 xTB，未填写能垒值。\n\n"
        "详细字段见 BARRIER_INTERFACE.md。几何检查失败的结构需先修复；失败计算记录不能填为零。\n",
        encoding="utf-8")
    shutil.copyfile(ROOT / "docs" / "BARRIER_INTERFACE.md", output / "BARRIER_INTERFACE.md")
    archive = shutil.make_archive(str(args.output / "external_evaluation"), "zip", output)
    progress("完成：Top3 与外部评估文件已导出")
    return {"output": str(output), "archive": archive}
