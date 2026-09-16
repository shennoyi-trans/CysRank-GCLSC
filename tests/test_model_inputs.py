import copy
import unittest
import numpy as np
import torch
from src import core, pipeline


class ModelInputTests(unittest.TestCase):
    def test_labels_do_not_enter_features_and_translation_invariance(self):
        torch.set_num_threads(2)
        row = pipeline.prepare(core.ROOT / 'data/raw', core.ROOT / 'data/provenance/parent_sequences.jsonl')[0]
        backbone = core.make_backbone(core.ROOT / 'models/pretrained/v_48_020.pt', torch.device('cpu'))
        clean = {k: copy.deepcopy(row[k]) for k in ('seq', 'xyz', *core.BIO)}
        original = core.features(backbone, clean, torch.device('cpu'))
        changed = copy.deepcopy(row)
        changed['ca_label'] = [1 - value for value in row['ca_label']]
        changed['supervision_weight'] = 500
        changed['experiment_outcome'] = 0
        self.assertTrue(torch.equal(original, core.features(backbone, changed, torch.device('cpu'))))
        clean['xyz'] = (np.asarray(clean['xyz']) + [10, -20, 30]).tolist()
        self.assertTrue(torch.allclose(original, core.features(backbone, clean, torch.device('cpu')), atol=1e-4, rtol=1e-4))

    def test_nonfinite_and_nested_features_rejected(self):
        row = core.load_jsonl(core.ROOT / 'data/examples/input.jsonl')[0]
        for invalid in (float('nan'), [0, 1, 0]):
            changed = copy.deepcopy(row)
            changed['sasa'][1] = invalid
            with self.assertRaises((ValueError, TypeError)):
                core.validate(changed)


if __name__ == '__main__':
    unittest.main()
