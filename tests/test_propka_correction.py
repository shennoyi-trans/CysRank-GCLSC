import math
import unittest
from tools.score_propka_case import correct


class TestFixedPkaCorrection(unittest.TestCase):
    def test_reference_preserves_score(self):
        self.assertAlmostEqual(correct(0.8, 8.0)[2], 0.8)

    def test_adjusts_odds_not_probability(self):
        old, new, score, alpha = correct(0.8, 9.0)
        self.assertAlmostEqual(new-old, -0.5)
        self.assertAlmostEqual(score/(1-score), 4*math.exp(-0.5))
        self.assertAlmostEqual(alpha, 1/(1+10**1.5))

    def test_invalid_inputs_fail(self):
        for score,pka in [(0,9),(1,9),(0.5,float('nan'))]:
            with self.assertRaises(ValueError):correct(score,pka)
