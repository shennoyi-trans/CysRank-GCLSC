import json
import socket
import threading
import urllib.request
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import app
from app import Jobs


class InputWindowTests(unittest.TestCase):
    def test_auto_port_and_browser_open(self):
        opened = threading.Event()
        servers = []
        browser_urls = []
        real_server = app.ThreadingHTTPServer

        def create_server(*args, **kwargs):
            server = real_server(*args, **kwargs)
            servers.append(server)
            return server

        def open_browser(url):
            browser_urls.append(url)
            opened.set()

        with socket.socket() as occupied:
            occupied.bind(('127.0.0.1', 0))
            occupied.listen()
            port = occupied.getsockname()[1]
            with patch('sys.argv', ['app.py', '--port', str(port), '--auto-port', '--open-browser']), \
                    patch.object(app, 'ThreadingHTTPServer', side_effect=create_server), \
                    patch.object(app.webbrowser, 'open', side_effect=open_browser):
                thread = threading.Thread(target=app.main, daemon=True)
                thread.start()
                try:
                    self.assertTrue(opened.wait(10))
                    self.assertNotEqual(servers[0].server_port, port)
                    with urllib.request.urlopen(browser_urls[0] + '/api/config', timeout=5) as response:
                        self.assertEqual(json.load(response)['max_length'], 100)
                    with urllib.request.urlopen(browser_urls[0], timeout=5) as response:
                        self.assertIn('环化智造', response.read().decode('utf-8'))
                finally:
                    for server in servers:
                        server.shutdown()
                    thread.join(10)
                self.assertFalse(thread.is_alive())

    def test_normalization_validation_and_concurrency(self):
        with tempfile.TemporaryDirectory() as d, patch("app.threading.Thread"):
            jobs = Jobs(SimpleNamespace(output=Path(d), max_length=10))
            for sequence in (None, "AX", "A", "A" * 11):
                with self.assertRaises(ValueError):
                    jobs.start(sequence)
            identity = jobs.start(" kt\ntks ")
            self.assertEqual(jobs.snapshot(identity)["sequence"], "KTTKS")
            with self.assertRaises(RuntimeError):
                jobs.start("AA")

    def test_failed_worker_releases_slot_and_keeps_error_log(self):
        with tempfile.TemporaryDirectory() as d, patch("app.threading.Thread"), patch("app.subprocess.run", return_value=SimpleNamespace(returncode=1)):
            jobs = Jobs(SimpleNamespace(output=Path(d), max_length=10, device="cpu", fold_model="test",
                                        chunk_size=32, local_files_only=True))
            identity = jobs.start("AA")
            jobs._run(identity)
            self.assertEqual(jobs.snapshot(identity)["status"], "failed")
            self.assertIsInstance(jobs.start("AA"), str)


if __name__ == "__main__":
    unittest.main()
