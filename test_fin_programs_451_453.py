#!/usr/bin/env python3
"""Regression tests for local FIN Programs P451--P453."""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import unittest

import numpy as np

import fin_programs_445_447 as p445
import fin_programs_448_450 as p448
import fin_programs_451_453 as current


ROOT = Path(__file__).resolve().parent


class Programs451453Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = json.loads(current.RESULTS_PATH.read_text(encoding="utf-8"))

    def test_p451_diagonal_theorem_and_coherent_counterexample(self) -> None:
        result = self.data["P451"]
        diagonal = [Fraction(value) for value in result["diagonal_half_distance_interval"]]
        coherent = [Fraction(value) for value in result["coherent_half_distance_interval"]]
        self.assertLess(diagonal[1] - diagonal[0], Fraction(1, 10**35))
        self.assertGreater(coherent[0], diagonal[1])
        self.assertGreater(float(coherent[0] - diagonal[1]), 5.6e-4)
        self.assertEqual(result["algebraic_stationarity"]["exact_root_residual"], "0")
        self.assertGreater(float(Fraction(result["coherent_psd_exact_minors"][0])), 0)
        self.assertGreater(float(Fraction(result["coherent_psd_exact_minors"][1])), 0)
        self.assertEqual(result["coherent_recursion_residual"], 0.0)
        self.assertIn("Open", result["status"])

        normalizer, _ = current.p451_coherent_normalizer()
        self.assertGreater(float(np.linalg.eigvalsh(normalizer)[0]), 0.057)
        self.assertAlmostEqual(float(np.trace(normalizer)), 1.0, places=14)
        self.assertTrue(current.P451_WITNESS_PATH.is_file())

    def test_p452_full_simplex_symmetrization_and_bound(self) -> None:
        result = self.data["P452"]
        condition = [Fraction(value) for value in result["certified_condition_interval"]]
        self.assertGreater(condition[0], 0)
        self.assertLess(result["full_simplex_gap"], 1e-3)
        self.assertGreater(result["improvement_over_P448_upper"], 0.0023)

        rng = np.random.default_rng(451453)
        for probabilities in rng.dirichlet(np.ones(4), size=24):
            symmetric = (probabilities + probabilities[::-1]) / 2
            original = current.p452_component_values(probabilities)
            improved = current.p452_component_values(symmetric)
            self.assertTrue(np.all(original <= improved + 2e-13))
            self.assertLessEqual(
                p445.p446_numeric_objective(probabilities),
                p445.p446_numeric_objective(symmetric) + 2e-13,
            )

    def test_p453_global_minimum_tv_uniqueness(self) -> None:
        result = self.data["P453"]
        self.assertEqual(result["forced_support_atoms"], 7)
        self.assertEqual(result["weight_signs"], [-1, 1, 1, 1, 1, -1, 1])
        self.assertGreater(float(Fraction(result["minimum_pairwise_node_separation"])), 0.058)
        self.assertGreater(float(Fraction(result["absolute_vandermonde_determinant_lower"])), 1.49e-9)
        negative = [Fraction(value) for value in result["negative_mass_interval"]]
        self.assertLess(negative[1] - negative[0], Fraction(1, 10**18))
        self.assertIn("unique", result["status"])
        self.assertIn("null-cycle gauge", result["removal_test"])

    def test_artifacts_and_epistemic_boundaries(self) -> None:
        for path in (
            current.RESULTS_PATH,
            current.SUMMARY_PATH,
            current.P451_PATH,
            current.P451_WITNESS_PATH,
            current.P452_PATH,
            current.P453_PATH,
            current.FIGURE_PATH,
        ):
            self.assertTrue(path.is_file(), path)
        metadata = self.data["metadata"]
        self.assertTrue(metadata["local_only"])
        self.assertFalse(metadata["external_physical_evidence"])
        self.assertFalse(metadata["selector_discharged"])
        self.assertFalse(metadata["dimensional_source_exported"])
        self.assertFalse(metadata["legacy_strict_bridge_complete"])
        self.assertFalse(metadata["legacy_role_transfer_started"])


if __name__ == "__main__":
    unittest.main()
