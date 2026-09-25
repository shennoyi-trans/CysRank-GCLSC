import json
import tempfile
import unittest
import contextlib
import io
from pathlib import Path
from types import SimpleNamespace
from src import design


class InsertionTests(unittest.TestCase):
    def test_true_insertion_and_mapping(self):
        for seq in ["AG", "ACA", "ACCA", "GPGAG"]:
            rows = design.insertion_candidates(seq)
            self.assertEqual(len(rows), len(seq) + 1)
            for gap, r in enumerate(rows):
                self.assertEqual(len(r["sequence"]), len(seq) + 1)
                self.assertEqual(r["sequence"][gap], "C")
                self.assertEqual(r["sequence"][:gap] + r["sequence"][gap + 1:], seq)
                self.assertIsNone(r["candidate_to_parent_1based"][gap])
                for j, original in enumerate(r["candidate_to_parent_1based"]):
                    if original is not None:
                        self.assertEqual(r["sequence"][j], seq[original - 1])

    def test_duplicate_sequences_keep_site_mapping(self):
        rows = design.insertion_candidates("ACA")
        duplicates = [r for r in rows if r["sequence"] == "ACCA"]
        self.assertEqual([r["inserted_position_1based"] for r in duplicates], [2, 3])

    def test_prepare_manifest_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            output = Path(d) / "design"
            rows = design.prepare("ACA", "fixture", output)
            manifest = json.loads((output / "structure_manifest.template.json").read_text())
            self.assertEqual(len(manifest), 3)
            self.assertTrue(all(r["structure_file"] is None for r in manifest))
            self.assertEqual(sum(len(r["candidate_ids"]) for r in manifest), len(rows))
            self.assertTrue(all("prediction_score" not in r for r in rows))
            with self.assertRaises(FileExistsError):
                design.prepare("ACA", "fixture", output)

    def test_bad_input(self):
        for seq in ("", "A", "AxA", "AXA", "A C", None):
            with self.assertRaises(ValueError):
                design.insertion_candidates(seq)
        for identity in ("../escape", "x/y", "", "has space"):
            with self.assertRaises(ValueError):
                design.insertion_candidates("AAA", identity)

    def test_ranking_deduplicates_sequences_and_retains_low_scores(self):
        rows = [{"candidate_id": "a", "sequence": "ACCA", "eligible": True, "prediction_score": 0.4},
                {"candidate_id": "b", "sequence": "ACCA", "eligible": True, "prediction_score": 0.3},
                {"candidate_id": "c", "sequence": "CACA", "eligible": False, "prediction_score": 0.99},
                {"candidate_id": "d", "sequence": "ACAC", "eligible": True, "prediction_score": 0.2}]
        top = design.rank_candidates(rows)
        self.assertEqual([r["candidate_id"] for r in top], ["a", "d", "c"])
        self.assertEqual(top[-1]["selection_status"], "fallback")
        self.assertIn("成功可能性较低", top[-1]["warning"])

    def test_fallback_excludes_missing_and_nonfinite_scores(self):
        rows = [{"candidate_id": str(i), "sequence": seq, "eligible": False,
                 "prediction_score": score, "exclusion_reasons": ["previous_residue_proline"]}
                for i, (seq, score) in enumerate([("ACA", None), ("AAC", float("nan")),
                                                 ("CAA", .7), ("CCA", .8), ("ACC", .6)])]
        top = design.rank_candidates(rows)
        self.assertEqual([r["sequence"] for r in top], ["CCA", "CAA", "ACC"])
        self.assertTrue(all("Pro" in r["warning"] for r in top))

    def test_insufficient_unique_sequences_are_not_duplicated(self):
        rows = [{**r, "eligible": False, "prediction_score": .5} for r in design.insertion_candidates("CC")]
        self.assertEqual(len(design.rank_candidates(rows)), 1)


def write_fixture(path, sequence):
    """Artificial extended coordinates for interface tests only, not a fold model."""
    from Bio.PDB import Structure, Model, Chain, Residue, Atom, PDBIO
    import numpy as np
    root = Structure.Structure("fixture")
    model = Model.Model(0)
    chain = Chain.Chain("A")
    root.add(model)
    model.add(chain)
    number = 1
    for i, aa in enumerate(sequence):
        residue = Residue.Residue((" ", i+1, " "), {"A": "ALA", "C": "CYS"}[aa], " ")
        chain.add(residue)
        coords = {"N": (0., 0., 0.), "CA": (1.45, 0.15, 0.), "C": (2.97, 0., 0.),
                  "O": (2.97, 1.23, 0.), "CB": (1.45, 0.15, 1.53)}
        if aa == "C":
            coords["SG"] = (1.45, 0.15, 3.34)
        for name, xyz in coords.items():
            residue.add(Atom.Atom(name, np.array(xyz)+[i*4.3, 0, 0], 0., 1., " ", f"{name:>4}", number, element=name[0]))
            number += 1
    writer = PDBIO()
    writer.set_structure(root)
    writer.save(str(path))


class StructureDesignTests(unittest.TestCase):
    def test_structure_features_and_bad_sequence(self):
        from src import structure
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "synthetic.pdb"
            write_fixture(path, "ACAA")
            sample, evidence, chain = structure.features(path, "ACAA", 2, "A")
            self.assertEqual(sample["feature_mask"], [0, 1, 0, 0])
            self.assertTrue(evidence["geometry_passed"])
            self.assertGreater(sample["sasa"][1], 0)
            self.assertEqual(sample["is_disulfide"][1], 0)
            with self.assertRaises(ValueError):
                structure.features(path, "AACA", 3, "A")
            with self.assertRaises(ValueError):
                structure.features(path, "ACAA", 2, "B")
            chain[2].detach_child("SG")
            structure.write_chain(chain, path)
            with self.assertRaises(ValueError):
                structure.features(path, "ACAA", 2, "A")

    def test_overlap_and_disulfide_evidence(self):
        from src import structure
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "synthetic.pdb"
            write_fixture(path, "CACA")
            chain, residues = structure.read_peptide(path, "CACA", "A")
            chain[3]["SG"].coord = chain[1]["SG"].coord + [2., 0., 0.]
            # Deliberately malformed synthetic candidate must fail geometry checks.
            structure.write_chain(chain, path)
            sample, evidence, _ = structure.features(path, "CACA", 3, "A")
            self.assertEqual(sample["is_disulfide"][2], 1)
            self.assertFalse(evidence["geometry_passed"])

    def test_real_test_pdb_existing_cys_is_not_insertion(self):
        from src import structure
        from src.paths import ROOT
        sample, evidence, _ = structure.features(ROOT / "data/structures/test.pdb", "AGCKNFFWKTFTSC", 3, "A")
        self.assertEqual(sample["is_disulfide"][2], 1)
        self.assertEqual(evidence["other_cys_SG_distances"][0]["position_1based"], 14)

    def test_full_design_and_barrier_request(self):
        from src import core
        with tempfile.TemporaryDirectory() as d, contextlib.redirect_stdout(io.StringIO()):
            root = Path(d)
            preparation = root / "prepare"
            design.prepare("AAAA", "synthetic_only", preparation)
            entries = json.loads((preparation / "structure_manifest.template.json").read_text())
            for i, entry in enumerate(entries):
                pdb = preparation / f"fixture{i}.pdb"
                write_fixture(pdb, entry["sequence"])
                entry.update(structure_file=pdb.name, chain_id="A", structure_method="synthetic unit-test fixture",
                             structure_notes="Not a real prediction or experimental structure")
            manifest = preparation / "structures.json"
            manifest.write_text(json.dumps(entries))
            args = SimpleNamespace(candidates=preparation / "candidates.jsonl", structures=manifest,
                                   checkpoint=core.ROOT / "models/site_predictor.pt", output=root / "result",
                                   top_k=3, allow_partial=False)
            design.score(args)
            summary = json.loads((args.output / "run.json").read_text())
            self.assertEqual(summary["candidate_count"], 5)
            self.assertEqual(summary["unique_top_sequences"], 3)
            rows = core.load_jsonl(args.output / "all_candidates.jsonl")
            self.assertIn("free_N_terminal_Cys_has_no_preceding_backbone_amide", rows[0]["exclusion_reasons"])
            requests = core.load_jsonl(args.output / "barrier_requests.jsonl")
            self.assertEqual(len(requests), 3)
            self.assertTrue(all(len(e["sample"]["seq"]) == 5 and e["barrier"]["value"] is None for e in requests))
            with self.assertRaises(FileExistsError):
                design.score(args)
            manifest.write_text(json.dumps(entries[:1]))
            args.output = root / "partial"
            with self.assertRaises(ValueError):
                design.score(args)
            args.allow_partial = True
            design.score(args)
            summary = json.loads((args.output / "run.json").read_text())
            self.assertEqual(summary["unique_top_sequences"], 1)
            self.assertEqual(summary["fallback_count"], 1)
            self.assertEqual(summary["scored_count"], 1)
            self.assertTrue((args.output / "barrier_requests.jsonl").exists())
            top = core.load_jsonl(args.output / "top3.jsonl")
            self.assertIn("成功可能性较低", top[0]["warning"])
            self.assertIn("成功可能性较低", (args.output / "report.md").read_text(encoding="utf-8"))
