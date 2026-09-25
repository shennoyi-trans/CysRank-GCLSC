import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from app import Jobs


class InputWindowTests(unittest.TestCase):
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
