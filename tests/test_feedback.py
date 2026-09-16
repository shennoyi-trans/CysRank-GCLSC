import copy
import contextlib
import io
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from src import pipeline as run


class FeedbackTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rows = run.prepare(run.ROOT / 'data/raw', run.ROOT / 'data/provenance/parent_sequences.jsonl')

    def event(self, **changes):
        sample = {k: copy.deepcopy(self.rows[0][k]) for k in ('seq', 'xyz', *run.core.BIO)}
        sample['seq'] = 'ACA'
        e = dict(observation_id='fixture-only-001', evidence_type='wet_experiment', split='train',
                 parent_group='fixture_parent', assay_id='fixture_assay', conditions={'test': 'synthetic unit-test fixture, not biological data'},
                 provenance='unit_test_only', sample=sample, position_1based=2, outcome=0,
                 reviewed=True, cyclization_efficiency_fraction=0.1)
        e.update(changes)
        return e

    def merge(self, events):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'feedback.jsonl'
            run.core.write_jsonl(path, events)
            return run.merge_feedback(self.rows, path)

    def test_main_only_masks(self):
        self.assertEqual(len(self.rows), 16)
        self.assertEqual(sum(sum(r['ca_label']) for r in self.rows), 5)
        self.assertEqual(sum(sum(r['confirmed_negative_mask']) for r in self.rows), 14)
        self.assertEqual(sum(sum(r['label_mask']) for r in self.rows), 19)
        self.assertTrue(all(r['supervision_weight'] == 1 for r in self.rows))
        self.assertTrue(all('低' not in r['source_file'] for r in self.rows))

    def test_rules_no_sasa_or_nearby_exclusion(self):
        row = copy.deepcopy(self.rows[0])
        for sasa in (0, 6.32, 150, 400):
            row['sasa'][1] = sasa
            row['nearby_atoms'][1] = 999
            self.assertEqual(run.reasons(row, 1, run.RULES), [])
        row['prev_is_pro'][1] = 1
        self.assertEqual(run.reasons(row, 1, run.RULES), ['previous_residue_proline'])
        row['is_disulfide'][1] = 1
        self.assertEqual(len(run.reasons(row, 1, run.RULES)), 2)

    def test_simulation_and_unreviewed_not_supervised(self):
        for event in (self.event(evidence_type='simulation'), self.event(reviewed=False), self.event(split='holdout')):
            rows, stats = self.merge([event])
            self.assertEqual(len(rows), 16)
            self.assertEqual(stats['accepted'], 0)

    def test_reviewed_wet_accepted_with_explicit_outcome(self):
        rows, stats = self.merge([self.event()])
        self.assertEqual(stats['accepted'], 1)
        self.assertEqual(rows[-1]['ca_label'], [0, 0, 0])
        self.assertEqual(rows[-1]['label_mask'], [0, 1, 0])
        self.assertEqual(rows[-1]['confirmed_negative_mask'], [0, 1, 0])

    def test_invalid_feedback_and_duplicates_rejected(self):
        for event in (self.event(outcome=None), self.event(position_1based=1),
                      self.event(cyclization_efficiency_fraction=75), self.event(conditions={}),
                      self.event(evidence_type='predicted')):
            with self.assertRaises(ValueError):
                run.validate_event(event)
        with self.assertRaises(ValueError):
            self.merge([self.event(), self.event()])
        with self.assertRaises(ValueError):
            self.merge([self.event(), self.event(observation_id='second')])

    def test_holdout_parent_and_sequence_leakage_rejected(self):
        with self.assertRaises(ValueError):
            self.merge([self.event(split='holdout', parent_group=self.rows[0]['parent_group'])])
        sample = {k: self.rows[0][k] for k in ('seq', 'xyz', *run.core.BIO)}
        with self.assertRaises(ValueError):
            self.merge([self.event(split='holdout', sample=sample)])

    def test_single_class_metrics_not_presented_as_accuracy(self):
        m = run.evaluate([1, 1], [0.8, 0.2])
        self.assertEqual(m['recall'], 0.5)
        self.assertIsNone(m['precision'])
        self.assertIsNone(m['accuracy'])

    def test_feedback_ingest_retrain_predict_end_to_end(self):
        # Synthetic fixtures live only in a temporary run; never in the delivered model.
        with tempfile.TemporaryDirectory() as directory, contextlib.redirect_stdout(io.StringIO()):
            folder = Path(directory)
            incoming = folder / 'incoming.jsonl'
            archive = folder / 'archive.jsonl'
            run.core.write_jsonl(incoming, [self.event(), self.event(
                observation_id='simulation-fixture', evidence_type='simulation')])
            run.ingest(SimpleNamespace(input=incoming, output=archive))
            output = folder / 'test_run'
            run.train(SimpleNamespace(source=run.ROOT / 'data/raw',
                parents=run.ROOT / 'data/provenance/parent_sequences.jsonl',
                pretrained=run.ROOT / 'models/pretrained/v_48_020.pt',
                feedback=archive, epochs=1, seed=1, output=output))
            summary = run.json.loads((output / 'summary.json').read_text())
            self.assertEqual(summary['training_sites'], 20)
            self.assertEqual(summary['feedback']['accepted'], 1)
            self.assertEqual(summary['feedback']['simulation'], 1)
            prediction = folder / 'predictions.jsonl'
            run.predict(SimpleNamespace(checkpoint=output / 'final_model/site_predictor.pt',
                input=output / 'prepared/training.jsonl', output=prediction))
            self.assertEqual(len(run.core.load_jsonl(prediction)), 20)


if __name__ == '__main__':
    unittest.main()
