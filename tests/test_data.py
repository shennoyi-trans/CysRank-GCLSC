"""Regression checks for the explicitly requested continued Pro exclusion."""
import unittest
from src import pipeline as run


class DataRevisionTests(unittest.TestCase):
    def test_exact_exclusion_and_backup(self):
        source = run.ROOT / 'data/raw'
        before = run.core.load_jsonl(run.ROOT / 'data/provenance/positive_before_exclusion.jsonl')
        after = run.core.load_jsonl(source / '训练正例.jsonl')
        self.assertEqual(len(before), 9)
        self.assertEqual(sum(r['seq'] == 'VNRRPCFSALT' for r in before), 1)
        self.assertEqual(after, [r for r in before if r['seq'] != 'VNRRPCFSALT'])

    def test_supervision_and_disclosed_test_overlap(self):
        rows = run.prepare(run.ROOT / 'data/raw',
                           run.ROOT / 'data/provenance/parent_sequences.jsonl')
        self.assertEqual(len(rows), 19)
        self.assertEqual(sum(sum(r['label_mask']) for r in rows), 22)
        self.assertEqual(sum(sum(r['ca_label']) for r in rows), 8)
        self.assertEqual(sum(sum(r['confirmed_negative_mask']) for r in rows), 14)
        self.assertEqual({r['seq'] for r in rows if r['seq'] in 'AGCKNFFWKTFTSC'},
                         {'AGCKNFFW', 'KTFTSC'})
        self.assertTrue(all(r['supervision_weight'] == 1 for r in rows))
        for r in rows:
            for i, y in enumerate(r['ca_label']):
                if y:
                    self.assertEqual(run.reasons(r, i, run.RULES), [])
        self.assertFalse(run.RULES['sasa_hard_filter'])

    def test_revision_hashes_and_unknown_labels(self):
        revision = run.json.loads((run.ROOT / 'data/provenance/revision_20260924.json').read_text(encoding='utf-8'))
        for path, digest in revision['active_file_sha256'].items():
            self.assertEqual(run.core.sha(run.ROOT / path), digest, path)
        external = run.core.load_jsonl(run.ROOT / 'data/external/validation.jsonl')
        self.assertEqual(len(external), 248)
        self.assertTrue(all(not any(r['ca_label']) for r in external))
        examples = run.core.load_jsonl(run.ROOT / 'data/examples/input.jsonl')
        self.assertEqual(sum(sum(r['feature_mask']) for r in examples), 22)
        self.assertTrue(all('ca_label' not in r for r in examples))


if __name__ == '__main__':
    unittest.main()
