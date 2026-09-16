"""Resolve every user-supplied relative path from the repository root."""
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]

def root_path(value):
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (ROOT / path).resolve()
