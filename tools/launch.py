"""Windows one-click setup; stdlib-only until dependencies have been installed."""
import argparse
import importlib.metadata
import json
from pathlib import Path
import platform
import subprocess
import sys
import uuid

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
REQUIREMENTS = ROOT / "requirements.txt"


def required_versions():
    return dict(line.strip().split("==", 1) for line in REQUIREMENTS.read_text(encoding="utf-8").splitlines()
                if line.strip() and not line.lstrip().startswith(("#", "-")))


def missing_requirements():
    missing = []
    for name, version in required_versions().items():
        try:
            installed = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            installed = None
        if installed != version:
            missing.append(f"{name}=={version}")
    return missing


def probe():
    compatible = (sys.version_info[:2] == (3, 13) and sys.maxsize > 2**32
                  and platform.python_implementation() == "CPython"
                  and not bool(__import__('sysconfig').get_config_var('Py_GIL_DISABLED')))
    missing = missing_requirements()
    # Prefer reusing the multi-GB Torch installation over many small packages.
    score = len(required_versions()) - len(missing)
    if not any(item.startswith("torch==") for item in missing):
        score += 100
    return {"compatible": compatible, "score": score, "missing": missing, "python": sys.executable}


def verify_assets():
    for relative in ("app.py", "design.py", "web/index.html", "docs/BARRIER_INTERFACE.md",
                     "models/site_predictor.pt", "models/pretrained/v_48_020.pt"):
        path = ROOT / relative
        if not path.is_file() or not path.stat().st_size:
            raise RuntimeError(f"项目文件缺失：{relative}。请完整解压项目后重试。")


def managed_python():
    """Reuse project environments; put additions over external packages in a local venv."""
    prefix = Path(sys.prefix).resolve()
    if sys.prefix != sys.base_prefix and prefix.is_relative_to(ROOT):
        return Path(sys.executable)
    folder = ROOT / ".runtime" / "venv-313"
    if folder.exists():
        # Never overwrite a copied/broken environment or an unrelated directory.
        folder = folder.with_name(f"venv-313-{uuid.uuid4().hex[:8]}")
    print(f"创建项目环境并复用当前 Python 已有的包：{folder}", flush=True)
    subprocess.run([sys.executable, "-m", "venv", "--system-site-packages", str(folder)], check=True)
    return folder / "Scripts" / "python.exe"


def install_dependencies():
    missing = missing_requirements()
    health = subprocess.run([sys.executable, "-m", "pip", "check"], capture_output=True, text=True)
    if missing or health.returncode:
        print("补齐依赖：" + (", ".join(missing) or "修复间接依赖"), flush=True)
        if health.returncode:
            print(health.stdout or health.stderr, flush=True)
        subprocess.run([sys.executable, "-m", "ensurepip", "--upgrade"], check=True)
        subprocess.run([sys.executable, "-m", "pip", "install", "--disable-pip-version-check",
                        "--prefer-binary", "-r", str(REQUIREMENTS)], check=True)
    else:
        print("完整依赖版本已满足，复用现有安装。", flush=True)
    subprocess.run([sys.executable, "-m", "pip", "check"], check=True)
    # Use a fresh process so metadata and DLLs are not cached from before installation.
    subprocess.run([sys.executable, str(Path(__file__).resolve()), "--verify-imports"], check=True)


def verify_imports():
    import numpy
    import torch
    from Bio.PDB import PDBParser
    from Bio.PDB.SASA import ShrakeRupley
    from transformers import EsmForProteinFolding
    import accelerate
    from src import core, structure
    # Exercise compiled Torch and NumPy modules rather than checking package names only.
    assert torch.from_numpy(numpy.ones(2)).sum().item() == 2
    print("页面预测、PDB 解析和评分依赖检查通过。", flush=True)


def prepare_model(offline=False):
    from huggingface_hub import hf_hub_download
    from huggingface_hub.errors import LocalEntryNotFoundError
    from src.folding import MODEL_ID, MODEL_REVISION
    names = ("config.json", "pytorch_model.bin")
    project_cache = ROOT / "models" / "folding-cache"
    # Also reuse a complete standard HF cache (including HF_HOME/HF_HUB_CACHE overrides).
    for cache in (str(project_cache), None):
        try:
            files = [Path(hf_hub_download(MODEL_ID, name, revision=MODEL_REVISION,
                                         cache_dir=cache, local_files_only=True)) for name in names]
            if all(file.is_file() and file.stat().st_size for file in files):
                print(f"复用 ESMFold 缓存：{files[0].parent}", flush=True)
                return files[0].parent
        except LocalEntryNotFoundError:
            pass
    if offline:
        raise RuntimeError("ESMFold 缓存尚不完整；正常双击 start.cmd 将自动下载。")
    print("下载 ESMFold 权重（约 8.44 GB）；已有完整文件会复用，中断后可重试。", flush=True)
    files = [Path(hf_hub_download(MODEL_ID, name, revision=MODEL_REVISION,
                                 cache_dir=str(project_cache))) for name in names]
    return files[0].parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probe", action="store_true")
    parser.add_argument("--verify-imports", action="store_true")
    parser.add_argument("--prepared-env", action="store_true")
    parser.add_argument("--check-only", action="store_true")
    parser.add_argument("--no-browser", action="store_true")
    parser.add_argument("--port", type=int, default=8765)
    args = parser.parse_args()
    if args.probe:
        print(json.dumps(probe()))
        return 0
    if args.verify_imports:
        verify_imports()
        return 0
    if not probe()["compatible"]:
        raise RuntimeError("需要 CPython 3.13 x64，请通过 start.cmd 自动配置。")
    if not 0 <= args.port <= 65535:
        raise ValueError("端口必须在 0–65535 之间。")
    verify_assets()
    if args.check_only:
        missing = missing_requirements()
        if missing:
            raise RuntimeError("需要安装：" + ", ".join(missing))
        subprocess.run([sys.executable, "-m", "pip", "check"], check=True)
        verify_imports()
        prepare_model(offline=True)
        print(f"检查通过；Python：{sys.executable}")
        return 0
    if not args.prepared_env:
        python = managed_python()
        command = [str(python), str(Path(__file__).resolve()), "--prepared-env", "--port", str(args.port)]
        if args.no_browser:
            command.append("--no-browser")
        return subprocess.call(command)
    install_dependencies()
    model = prepare_model()
    print("环境和模型已就绪，启动后端。关闭此窗口或按 Ctrl+C 可停止服务。", flush=True)
    command = [sys.executable, "-u", str(ROOT / "app.py"), "--port", str(args.port),
               "--fold-model", str(model), "--local-files-only", "--auto-port"]
    if not args.no_browser:
        command.append("--open-browser")
    return subprocess.call(command, cwd=ROOT)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("\n服务已停止。")
        raise SystemExit(0)
    except Exception as exc:
        print(f"\n启动失败：{exc}\n修复网络或文件问题后重新双击 start.cmd；已安装的依赖和下载缓存会保留。", file=sys.stderr)
        raise SystemExit(1)
