"""Local sequence-input window: python app.py, then open http://127.0.0.1:8765."""
import argparse
import json
import mimetypes
import subprocess
import sys
import threading
import uuid
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlparse

from src.design import insertion_candidates
from src.paths import ROOT, root_path


class Jobs:
    def __init__(self, options):
        self.options = options
        self.rows = {}
        self.lock = threading.Lock()

    def start(self, sequence):
        if not isinstance(sequence, str):
            raise ValueError("请输入多肽序列")
        sequence = "".join(sequence.split()).upper()
        insertion_candidates(sequence)
        if len(sequence) > self.options.max_length:
            raise ValueError(f"当前配置最多支持 {self.options.max_length} 个残基")
        with self.lock:
            if any(row["status"] == "running" for row in self.rows.values()):
                raise RuntimeError("已有预测正在运行，请等待完成后再提交")
            identity = uuid.uuid4().hex
            folder = self.options.output / identity
            folder.mkdir(parents=True)
            self.rows[identity] = {"id": identity, "sequence": sequence, "status": "running", "folder": folder}
        threading.Thread(target=self._run, args=(identity,), daemon=True).start()
        return identity

    def _run(self, identity):
        row = self.rows[identity]
        command = [sys.executable, "-u", str(ROOT / "design.py"), "run", "--sequence", row["sequence"],
                   "--record-id", "peptide", "--output", str(row["folder"] / "design"),
                   "--device", self.options.device, "--fold-model", self.options.fold_model,
                   "--chunk-size", str(self.options.chunk_size)]
        if self.options.local_files_only:
            command.append("--local-files-only")
        try:
            import os
            env = {**os.environ, "PYTHONIOENCODING": "utf-8"}
            with (row["folder"] / "progress.log").open("w", encoding="utf-8") as log:
                result = subprocess.run(command, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT,
                                        creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
            with self.lock:
                row["status"] = "completed" if result.returncode == 0 else "failed"
        except Exception as error:
            with self.lock:
                row.update(status="failed", error=str(error))

    def snapshot(self, identity):
        with self.lock:
            row = dict(self.rows[identity])
        folder = row.pop("folder")
        log = folder / "progress.log"
        if log.exists():
            with log.open("rb") as stream:
                stream.seek(max(0, log.stat().st_size - 12000))
                row["log"] = stream.read().decode("utf-8", errors="replace")
        if row["status"] == "completed":
            result = folder / "design" / "results"
            row["top3"] = [json.loads(line) for line in (result / "top3.jsonl").read_text(encoding="utf-8").splitlines() if line]
            row["summary"] = json.loads((result / "run.json").read_text(encoding="utf-8"))
        return row


def handler_for(jobs):
    class Handler(BaseHTTPRequestHandler):
        def send(self, data, status=200, content_type="application/json; charset=utf-8", attachment=None):
            if not isinstance(data, bytes):
                data = json.dumps(data, ensure_ascii=False).encode("utf-8")
            self.send_response(status)
            self.send_header("Content-Type", content_type)
            self.send_header("Content-Length", str(len(data)))
            self.send_header("Cache-Control", "no-store")
            self.send_header("X-Content-Type-Options", "nosniff")
            if attachment:
                self.send_header("Content-Disposition", f'attachment; filename="{attachment}"')
            self.end_headers()
            self.wfile.write(data)

        def local_request(self):
            allowed = {f"127.0.0.1:{self.server.server_port}", f"localhost:{self.server.server_port}"}
            if self.headers.get("Host") not in allowed:
                self.send({"error": "Local access only"}, 403)
                return False
            origin = self.headers.get("Origin")
            if origin and origin not in {"http://" + value for value in allowed}:
                self.send({"error": "Same-origin requests only"}, 403)
                return False
            return True

        def do_GET(self):
            if not self.local_request():
                return
            path = urlparse(self.path).path
            if path == "/":
                return self.send((ROOT / "web" / "index.html").read_bytes(), content_type="text/html; charset=utf-8")
            if path == "/api/config":
                return self.send({"max_length": jobs.options.max_length, "device": jobs.options.device})
            parts = path.strip("/").split("/")
            if len(parts) == 3 and parts[:2] == ["api", "jobs"]:
                try:
                    return self.send(jobs.snapshot(parts[2]))
                except KeyError:
                    return self.send({"error": "任务不存在，请重新提交"}, 404)
            if len(parts) >= 3 and parts[0] == "downloads":
                row = jobs.rows.get(parts[1])
                if row and row["status"] == "completed":
                    base = (row["folder"] / "design").resolve()
                    file = (base / "/".join(parts[2:])).resolve()
                    if file.is_relative_to(base) and file.is_file() and file.suffix in (".zip", ".pdb", ".csv", ".jsonl", ".md"):
                        return self.send(file.read_bytes(), content_type=mimetypes.guess_type(file.name)[0] or "application/octet-stream",
                                         attachment=file.name)
            self.send({"error": "Not found"}, 404)

        def do_POST(self):
            if not self.local_request():
                return
            if self.path != "/api/jobs":
                return self.send({"error": "Not found"}, 404)
            try:
                size = int(self.headers.get("Content-Length", 0))
                if not 0 < size <= 16384:
                    raise ValueError("输入内容过长或为空")
                payload = json.loads(self.rfile.read(size))
                if not isinstance(payload, dict):
                    raise ValueError("需要 JSON 对象")
                identity = jobs.start(payload.get("sequence"))
                self.send({"id": identity}, 202)
            except (ValueError, TypeError) as error:
                self.send({"error": str(error)}, 400)
            except RuntimeError as error:
                self.send({"error": str(error)}, 409)
            except OSError as error:
                self.send({"error": f"无法创建结果目录，请检查运行权限：{error}"}, 500)

    return Handler


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--output", type=root_path, default=ROOT / "results" / "web_jobs")
    parser.add_argument("--device", choices=["auto", "cpu", "cuda"], default="auto")
    parser.add_argument("--fold-model", default="facebook/esmfold_v1")
    parser.add_argument("--chunk-size", type=int, default=32)
    parser.add_argument("--max-length", type=int, default=100)
    parser.add_argument("--local-files-only", action="store_true")
    options = parser.parse_args()
    if options.max_length < 2 or options.chunk_size < 1:
        parser.error("max-length must be >= 2 and chunk-size must be positive")
    server = ThreadingHTTPServer(("127.0.0.1", options.port), handler_for(Jobs(options)))
    print(f"CysRank input window: http://127.0.0.1:{server.server_port}", flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
