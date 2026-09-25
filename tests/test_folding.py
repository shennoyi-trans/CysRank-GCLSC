import contextlib
import io
import json
import tempfile
import unittest
import zipfile
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from src import design, folding
from src.paths import ROOT
from test_design import write_fixture


class FixturePredictor:
    """TEST ONLY. Never selected by the production application."""
    metadata = {"method": "synthetic test fixture, not a prediction"}

    def __init__(self, *args, **kwargs):
        self.calls = []

    def predict(self, sequence):
        self.calls.append(sequence)
        with tempfile.TemporaryDirectory() as d:
            file = Path(d) / "fixture.pdb"
            write_fixture(file, sequence)
            return file.read_text()


class FoldingTests(unittest.TestCase):
    def test_predicts_unique_sequences_and_rejects_mismatch(self):
        with tempfile.TemporaryDirectory() as d:
            preparation = Path(d) / "prepare"
            design.prepare("ACA", "test", preparation)
            predictor = FixturePredictor()
            manifest = folding.predict_manifest(preparation, predictor, progress=lambda _: None)
            self.assertEqual(len(predictor.calls), 3)
            self.assertEqual(len(json.loads(manifest.read_text())), 3)
            other = Path(d) / "bad"
            design.prepare("AA", "test", other)
            predictor.predict = lambda _: "ATOM fake"
            with self.assertRaises(ValueError):
                folding.predict_manifest(other, predictor, progress=lambda _: None)

    def test_orchestration_exports_portable_package_and_fallback(self):
        from src.barrier import sample_digest
        with tempfile.TemporaryDirectory() as d, patch.object(folding, "ESMFoldPredictor", FixturePredictor), contextlib.redirect_stdout(io.StringIO()):
            args = SimpleNamespace(output=Path(d) / "run", sequence="AA", record_id="fixture",
                                   fold_model="test", device="cpu", chunk_size=32, local_files_only=True,
                                   checkpoint=ROOT / "models/site_predictor.pt", top_k=3)
            result = folding.run_design(args, progress=lambda _: None)
            output = Path(result["output"])
            top = [json.loads(line) for line in (output / "top3.jsonl").read_text(encoding="utf-8").splitlines()]
            self.assertEqual(len(top), 3)
            self.assertEqual(top[-1]["selection_status"], "fallback")
            requests = [json.loads(line) for line in (output / "external_barrier_requests.jsonl").read_text(encoding="utf-8").splitlines()]
            for row in requests:
                self.assertEqual(row["barrier"]["sample_sha256"], sample_digest(row["sample"]))
                self.assertTrue((output / row["barrier"]["structure_reference"]).is_file())
                self.assertIsNone(row["barrier"]["value"])
            with zipfile.ZipFile(result["archive"]) as archive:
                self.assertIn("top3.csv", archive.namelist())
                self.assertIn("EXTERNAL_EVALUATION.md", archive.namelist())
            with self.assertRaises(FileExistsError):
                folding.run_design(args)


if __name__ == "__main__":
    unittest.main()
