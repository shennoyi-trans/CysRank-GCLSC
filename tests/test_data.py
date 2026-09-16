"""Regression checks for the explicitly requested ZQY8 exclusion."""
import unittest
from src import pipeline as run


class ZQY8Tests(unittest.TestCase):
    def test_exact_exclusion_and_backup(self):
        source = run.ROOT / 'data/raw'
        before = run.core.load_jsonl(run.ROOT / 'data/provenance/positive_before_exclusion.jsonl')
        after = run.core.load_jsonl(source / '训练正例.jsonl')
        self.assertEqual(sum(r['seq'] == 'VNRRPCFSALT' for r in before), 1)
        self.assertEqual(after, [r for r in before if r['seq'] != 'VNRRPCFSALT'])

    def test_supervision_and_test_separation(self):
        rows = run.prepare(run.ROOT / 'data/raw',
                           run.ROOT / 'data/provenance/parent_sequences.jsonl')
        self.assertEqual(len(rows), 16)
        self.assertEqual(sum(sum(r['label_mask']) for r in rows), 19)
        self.assertEqual(sum(sum(r['ca_label']) for r in rows), 5)
        self.assertEqual(sum(sum(r['confirmed_negative_mask']) for r in rows), 14)
        self.assertFalse(any(r['seq'] in 'AGCKNFFWKTFTSC' or
                             'AGCKNFFWKTFTSC' in r['seq'] for r in rows))
        self.assertTrue(all(r['supervision_weight'] == 1 for r in rows))
        for r in rows:
            for i, y in enumerate(r['ca_label']):
                if y:
                    self.assertEqual(run.reasons(r, i, run.RULES), [])
        self.assertFalse(run.RULES['sasa_hard_filter'])


if __name__ == '__main__':
    unittest.main()
