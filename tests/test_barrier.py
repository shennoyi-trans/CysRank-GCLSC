"""Synthetic fixtures only: no calculated or experimental claims."""
import contextlib
import copy
import io
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from src import barrier, core, pipeline


class BarrierTests(unittest.TestCase):
    def event(self, identity="a", seq="ACA", split="train", value=10.0):
        sample = {"seq": seq, "xyz": [[0, 0, 0], [3.8, 0, 0], [7, 2, 0]],
                  "is_disulfide": [0, 0, 0], "sasa": [0, 80, 0],
                  "prev_is_pro": [0, 0, 0], "nearby_atoms": [0, 10, 0]}
        return {"observation_id": identity, "evidence_type": "simulation", "split": split,
                "parent_group": identity, "assay_id": "synthetic-unit-test", "conditions": {"fixture": True},
                "provenance": "synthetic test only, not a real calculation", "sample": sample,
                "position_1based": 2, "outcome": None, "reviewed": True,
                "barrier": {"schema_version": barrier.VERSION, "status": "completed", "value": value,
                            "unit": "kcal/mol", "quantity": "delta_g_dagger", "quality_passed": True,
                            "quality_notes": "synthetic test", "structure_reference": "synthetic.xyz",
                            "sample_sha256": barrier.sample_digest(sample), "charge": 0, "multiplicity": 1,
                            "protocol": {**{k: "synthetic" for k in barrier.PROTOCOL_FIELDS},
                                         "temperature_kelvin": 298.15}}}

    def test_units_signed_values_and_invalid_numbers(self):
        e = self.event()
        for unit, factor in barrier.FACTORS.items():
            e["barrier"].update(unit=unit, value=10 / factor)
            self.assertAlmostEqual(barrier.validate_barrier(e), 10)
        e["barrier"].update(unit="kcal/mol", value=-1)
        self.assertEqual(barrier.validate_barrier(e), -1)
        for value in (True, float("nan"), float("inf"), "10", None):
            e["barrier"]["value"] = value
            with self.assertRaises(ValueError):
                pipeline.validate_event(e)

    def test_bad_contracts_and_sample_tampering_rejected(self):
        for key, value in (("unit", "kcal"), ("quantity", "reaction_energy"),
                           ("status", "pending"), ("quality_passed", "true"),
                           ("charge", 0.5), ("multiplicity", 0), ("sample_sha256", "wrong")):
            e = self.event()
            e["barrier"][key] = value
            with self.assertRaises(ValueError):
                pipeline.validate_event(e)
        e = self.event()
        e["sample"]["sasa"][1] += 1
        with self.assertRaises(ValueError):
            pipeline.validate_event(e)
        for change in ({"evidence_type": "wet_experiment"}, {"outcome": 1}):
            e = self.event()
            e.update(change)
            with self.assertRaises(ValueError):
                pipeline.validate_event(e)

    def test_failed_unreviewed_and_quality_failed_excluded(self):
        for mode in ("failed", "unreviewed", "quality"):
            e = self.event()
            if mode == "failed":
                e["barrier"].update(status="failed", value=None, quality_passed=False, failure_reason="SCF failed")
            elif mode == "unreviewed":
                e["reviewed"] = False
            else:
                e["barrier"]["quality_passed"] = False
            rows, stats = barrier.collect([e])
            self.assertEqual(rows, [])
            self.assertEqual(stats["ineligible"], 1)
        e["barrier"].update(status="failed", value=0, quality_passed=False, failure_reason="SCF failed")
        with self.assertRaises(ValueError):
            pipeline.validate_event(e)

    def test_duplicates_protocol_mixing_and_leakage(self):
        a = self.event()
        cases = [[a, copy.deepcopy(a)], [a, self.event("b")],
                 [a, self.event("b", split="holdout")]]
        b = self.event("b", "GCG", "holdout")
        b["parent_group"] = "a"
        cases.append([a, b])
        b = self.event("b", "GCG")
        b["barrier"]["quantity"] = "delta_e_dagger"
        cases.append([a, b])
        b = self.event("b", "GCG")
        b["barrier"]["protocol"]["method"] = "different"
        cases.append([a, b])
        for events in cases:
            with self.assertRaises(ValueError):
                barrier.collect(events)

    def test_scaling_and_intercept(self):
        x = np.array([[1., 2.], [3., 2.]])
        y = np.array([10., 20.])
        model = barrier.fit_ridge(x, y, 1.)
        np.testing.assert_allclose(model["mean"], [2., 2.])
        self.assertEqual(model["intercept"], 15.)
        self.assertTrue(np.isfinite(barrier.apply_regressor(model, x)).all())
        for alpha in (0, -1, float("nan")):
            with self.assertRaises(ValueError):
                barrier.fit_ridge(x, y, alpha)

    def test_end_to_end_archive_train_reload_predict_and_wet_isolation(self):
        with tempfile.TemporaryDirectory() as folder, contextlib.redirect_stdout(io.StringIO()):
            root = Path(folder)
            events = [self.event(), self.event("b", "GCG", value=20),
                      self.event("c", "SCS", "holdout", value=100)]
            incoming, archive = root / "incoming.jsonl", root / "archive.jsonl"
            core.write_jsonl(incoming, events)
            pipeline.ingest(SimpleNamespace(input=incoming, output=archive))
            wet_rows, stats = pipeline.merge_feedback([], archive)
            self.assertEqual(wet_rows, [])
            self.assertEqual(stats["accepted"], 0)
            pretrained = core.ROOT / "models/pretrained/v_48_020.pt"
            output = root / "trained"
            barrier.train(SimpleNamespace(input=archive, output=output, alpha=1., pretrained=pretrained))
            model = json.loads((output / "barrier_model.json").read_text())
            self.assertEqual(model["training_observations"], ["a", "b"])
            self.assertEqual(model["head"]["intercept"], 15.)
            # Changing a held-out target cannot change learned weights/scaling.
            events[-1]["barrier"]["value"] = 9000
            core.write_jsonl(incoming, events)
            barrier.train(SimpleNamespace(input=incoming, output=root / "second", alpha=1., pretrained=pretrained))
            second = json.loads((root / "second/barrier_model.json").read_text())
            self.assertEqual(model["head"], second["head"])
            sample = {**events[0]["sample"], "record_id": "new", "feature_mask": [0, 1, 0]}
            inputs = root / "predict.jsonl"
            core.write_jsonl(inputs, [sample])
            args = SimpleNamespace(input=inputs, output=root / "prediction.jsonl", model=output / "barrier_model.json",
                                   protocol=output / "protocol.json", pretrained=pretrained)
            barrier.predict(args)
            prediction = core.load_jsonl(args.output)[0]
            self.assertIn("predicted_barrier_kcal_mol", prediction)
            self.assertNotIn("probability", prediction)
            with self.assertRaises(FileExistsError):
                barrier.predict(args)
            bad = json.loads(args.protocol.read_text())
            bad["quantity"] = "delta_e_dagger"
            core.write_json(root / "bad.json", bad)
            args.protocol, args.output = root / "bad.json", root / "bad_prediction.jsonl"
            with self.assertRaises(ValueError):
                barrier.predict(args)
            barrier.requests(SimpleNamespace(input=inputs, output=root / "requests.jsonl"))
            request = core.load_jsonl(root / "requests.jsonl")[0]
            self.assertIsNone(request["barrier"]["value"])
            self.assertEqual(request["barrier"]["status"], "pending")
            self.assertEqual(request["barrier"]["sample_sha256"], barrier.sample_digest(sample))
            with self.assertRaises(ValueError):
                pipeline.validate_event(request)


if __name__ == "__main__":
    unittest.main()
