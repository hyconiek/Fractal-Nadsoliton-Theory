#!/usr/bin/env python3
"""Regression tests for local FIN Programs P448--P450."""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import unittest

import numpy as np

import fin_programs_435_436_440 as p435
import fin_programs_445_447 as p445
import fin_programs_448_450 as current


ROOT = Path(__file__).resolve().parent


class Programs448450Tests(unittest.TestCase):
    def test_p448_majorant_encloses_and_symmetrizes(self) -> None:
        rng = np.random.default_rng(448450)
        for probabilities in rng.dirichlet(np.ones(4), size=16):
            point = [(float(value), float(value)) for value in probabilities]
            interval = current.p448_majorant_float_interval(point)
            majorant = current.p448_majorant_numeric(probabilities)
            original = p445.p446_numeric_objective(probabilities)
            symmetric = (probabilities + probabilities[::-1]) / 2
            symmetric_majorant = current.p448_majorant_numeric(symmetric)
            self.assertLessEqual(interval[0], majorant)
            self.assertGreaterEqual(interval[1], majorant)
            self.assertLessEqual(original, majorant + 2e-13)
            self.assertLessEqual(majorant, symmetric_majorant + 2e-13)

    def test_p448_global_certificate_scope(self) -> None:
        data = json.loads((ROOT / current.RESULTS_PATH.name).read_text(encoding="utf-8"))
        result = data["P448"]
        self.assertLess(result["original_certified_lower"], result["original_full_simplex_global_upper"])
        self.assertLess(result["original_full_simplex_gap"], 0.00336)
        self.assertGreater(result["original_full_simplex_gap"], 0)
        self.assertIn("[Open] exact full-simplex globality", result["status"])
        self.assertIn("strictly larger", result["proof_boundary"])

    def test_p449_recursion_dimensions_and_echo_certificate(self) -> None:
        self.assertEqual([current.support_affine_dimension(n) for n in range(1, 6)], [1, 5, 21, 85, 341])
        echo_interval, ghz_interval = current.p449_echo_witness_interval()
        self.assertGreater(echo_interval[0], ghz_interval[1])
        self.assertGreater(float(echo_interval[0] - ghz_interval[1]), 0.0176)

        echo, xi3, tester_sum = current.p449_full_echo_extension()
        self.assertGreaterEqual(float(np.linalg.eigvalsh(xi3)[0]), 0.0249999999999998)
        xi3_trace = p435.partial_trace_square(xi3, (2, 2, 2, 2, 2), 0)
        self.assertLess(np.linalg.norm(xi3_trace - np.eye(16) / 4), 1e-14)
        support = [
            np.ravel_multi_index((a, a, b, b, c, c), (2, 2, 2, 2, 2, 2))
            for a in (0, 1) for b in (0, 1) for c in (0, 1)
        ]
        self.assertLess(
            np.linalg.norm(tester_sum[np.ix_(support, support)] - echo), 1e-14
        )
        self.assertTrue((ROOT / current.P449_WITNESS_PATH.name).is_file())

    def test_p450_exact_null_cycle_and_representation_obstruction(self) -> None:
        nodes, weights = current.null_cycle()
        for order in range(12):
            self.assertEqual(sum(w * x**order for x, w in zip(nodes, weights)), 0)
        data = json.loads((ROOT / current.RESULTS_PATH.name).read_text(encoding="utf-8"))
        result = data["P450"]
        self.assertTrue(result["same_sign_identical_node_split_invariant"])
        self.assertEqual(Fraction(result["maximum_exact_null_moment_residual"]), 0)
        self.assertGreater(result["largest_tested_risk_ratio"], 300)
        self.assertIn("unbounded", result["unboundedness_theorem"])

    def test_artifacts_and_epistemic_boundaries(self) -> None:
        for path in (
            current.RESULTS_PATH,
            current.SUMMARY_PATH,
            current.P448_PATH,
            current.P449_PATH,
            current.P449_WITNESS_PATH,
            current.P450_PATH,
            current.FIGURE_PATH,
        ):
            self.assertTrue(path.is_file(), path)
        data = json.loads(current.RESULTS_PATH.read_text(encoding="utf-8"))
        self.assertTrue(data["metadata"]["local_only"])
        self.assertFalse(data["metadata"]["external_physical_evidence"])
        self.assertFalse(data["metadata"]["selector_discharged"])
        self.assertFalse(data["metadata"]["dimensional_source_exported"])
        self.assertFalse(data["metadata"]["legacy_strict_bridge_complete"])


if __name__ == "__main__":
    unittest.main()
