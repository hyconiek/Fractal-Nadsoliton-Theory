#!/usr/bin/env python3
"""Regression and certificate tests for FIN Programs P435/P436/P440."""

from __future__ import annotations

from fractions import Fraction
import json
import math
from pathlib import Path
import unittest

import csv
import mpmath as mp
import numpy as np

import fin_programs_435_436_440 as batch


ROOT = Path(__file__).resolve().parent


class Programs435436440Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.results = json.loads(batch.RESULTS_PATH.read_text(encoding="utf-8"))

    def test_p435_exact_one_slot_primal_dual(self) -> None:
        result, _ = batch.program_435()
        expected = (1 + (4 / 5) * math.sin(math.pi / 8)) / 2
        self.assertAlmostEqual(result["one_slot_optimal_success"], expected, places=14)
        checks = result["one_slot_checks"]
        self.assertLess(abs(checks["primal_dual_gap"]), 1e-13)
        self.assertLess(checks["tester_normalization_defect"], 1e-13)
        self.assertLess(checks["dual_reduction_scalar_defect"], 1e-13)
        self.assertLess(result["two_slot_checks"]["comb_plus_causality_defect"], 1e-13)
        self.assertEqual(result["two_slot_checks"]["primal_real_variables_before_equalities"], 580)
        self.assertTrue(all(not installed for installed in result["solver_inventory"].values()))

    def test_p436_machine_readable_certificate_is_strict(self) -> None:
        result = self.results["P436"]
        candidate_lower = float(result["candidate_total_interval"][0])
        product_upper = float(result["product_total_interval"][1])
        ghz_upper = float(result["ghz_total_interval"][1])
        self.assertGreater(candidate_lower, product_upper)
        self.assertGreater(candidate_lower, ghz_upper)
        self.assertGreater(float(result["certified_gain_lower"]), 0.0225)
        self.assertEqual(sum(Fraction(value) for value in result["candidate_probabilities"]), 1)

    def test_p436_analytic_constant_enclosures(self) -> None:
        self.assertLess(batch.PI_INTERVAL[0], batch.PI_INTERVAL[1])
        self.assertLess(float(batch.PI_INTERVAL[1] - batch.PI_INTERVAL[0]), 1e-80)
        self.assertEqual(float(batch.PI_INTERVAL[0]), math.pi)
        self.assertEqual(float(batch.PI_INTERVAL[1]), math.pi)
        mp.mp.dps = 100
        for difference in range(-3, 4):
            low, high = batch.p436_sine_interval(difference)
            target = mp.sin(2 * mp.pi * difference / 15)
            low_mp = mp.mpf(low.numerator) / low.denominator
            high_mp = mp.mpf(high.numerator) / high.denominator
            self.assertLessEqual(low_mp, target)
            self.assertGreaterEqual(high_mp, target)

    def test_p440_probabilities_and_improvements(self) -> None:
        result = self.results["P440"]
        with batch.P440_PATH.open(encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        minimax = np.array([float(row["P440_minimax_probability"]) for row in rows])
        self.assertAlmostEqual(float(minimax.sum()), 1.0, places=13)
        self.assertTrue(np.all(minimax > 0))
        coefficients = result["worst_mse_coefficients"]
        self.assertLess(coefficients["P440_detector_minimax"], coefficients["P422_spectral_radius"])
        self.assertLess(coefficients["P422_spectral_radius"], coefficients["absolute_weight_baseline"])
        self.assertGreater(result["confidence_ledger"]["hoeffding_sufficient_attempts"], 1_000_000)

    def test_artifacts_and_local_only_boundary(self) -> None:
        for path in (
            batch.RESULTS_PATH, batch.SUMMARY_PATH, batch.P435_PATH,
            batch.P436_PATH, batch.P440_PATH, batch.COMB_INSTANCE_PATH,
            batch.FIGURE_PATH,
        ):
            self.assertTrue(path.is_file(), path)
            self.assertGreater(path.stat().st_size, 0)
        metadata = self.results["metadata"]
        self.assertTrue(metadata["local_only"])
        self.assertFalse(metadata["external_physical_evidence"])
        self.assertFalse(metadata["selector_discharged"])
        self.assertFalse(metadata["dimensional_source_exported"])
        self.assertFalse(metadata["legacy_strict_bridge_complete"])


if __name__ == "__main__":
    unittest.main()
