#!/usr/bin/env python3
"""Regression tests for FIN local Programs P445--P447."""

from __future__ import annotations

from fractions import Fraction
import json
import math
from pathlib import Path
import unittest

import numpy as np

import fin_programs_411_427 as old
import fin_programs_435_436_440 as prior
import fin_programs_445_447 as batch


ROOT = Path(__file__).resolve().parent


class Programs445447Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.results = json.loads(batch.RESULTS_PATH.read_text(encoding="utf-8"))

    def test_p445_exact_value_and_condition(self) -> None:
        result = self.results["P445"]
        expected = 8 * math.sqrt(2) / 25
        self.assertAlmostEqual(result["optimal_half_distance"], expected, places=15)
        self.assertAlmostEqual(result["optimal_success"], (1 + expected) / 2, places=15)
        self.assertGreater(result["condition_margin_float"], 0.06)
        self.assertLess(result["maximum_random_formula_defect"], 1e-13)
        self.assertEqual(result["optimal_history_law"], ["1/2", "0", "0", "1/2"])

    def test_p445_general_normalizer_formula(self) -> None:
        probabilities = np.array([0.2, 0.3, 0.3, 0.2])
        coherence = 0.1 + 0.05j
        self.assertLess(abs(coherence) ** 2, min(probabilities[0] * probabilities[1], probabilities[2] * probabilities[3]))
        q = 0.8; theta = math.pi / 8
        j_plus = prior.qubit_phase_dephasing_choi(q, theta, +1)
        j_minus = prior.qubit_phase_dephasing_choi(q, theta, -1)
        delta = (np.kron(j_plus, j_plus) - np.kron(j_minus, j_minus))[np.ix_([0, 3, 12, 15], [0, 3, 12, 15])]
        normalizer = np.zeros((4, 4), dtype=complex)
        normalizer[:2, :2] = [[probabilities[0], coherence], [coherence.conjugate(), probabilities[1]]]
        normalizer[2:, 2:] = [[probabilities[2], -coherence], [-coherence.conjugate(), probabilities[3]]]
        values, vectors = np.linalg.eigh(normalizer)
        root = (vectors * np.sqrt(values)) @ vectors.conj().T
        direct = 0.5 * np.sum(np.abs(np.linalg.eigvalsh(root @ delta @ root)))
        formula = batch.p445_reduced_distance(probabilities, q, theta, coherence)
        self.assertAlmostEqual(direct, formula, places=14)

    def test_p446_fast_block_formula_matches_original(self) -> None:
        rng = np.random.default_rng(4460)
        for probabilities in rng.dirichlet(np.ones(4), size=12):
            fast = batch.p446_numeric_objective(probabilities)
            direct = old.erasure_aware_distance(3, 2 * math.pi / 15, 0.8, 0.8, probabilities)
            self.assertAlmostEqual(fast, direct, places=13)

    def test_p446_scope_is_not_promoted(self) -> None:
        result = self.results["P446"]
        certificate = result["palindromic_certificate"]
        self.assertLess(certificate["optimality_gap"], 0.001)
        self.assertGreaterEqual(certificate["global_upper"], certificate["incumbent_lower"])
        self.assertEqual(result["palindromic_family"], "p(a)=(a,1/2-a,1/2-a,a), 0<=a<=1/2")
        self.assertIn("full symmetric Hamming-sector simplex remains open", result["boundary"])
        self.assertLess(result["full_simplex_audit"]["maximum_sampled_hessian_eigenvalue"], 0)
        self.assertLess(result["full_simplex_audit"]["maximum_best_reversal_defect"], 1e-7)

    def test_p447_complete_box_propagation(self) -> None:
        result = self.results["P447"]
        self.assertTrue(result["strict_order_certified"])
        self.assertLess(float(result["maximum_probability_width"]), 4e-20)
        baseline = [float(value) for value in result["mse_absolute_weight_interval"]]
        p422 = [float(value) for value in result["mse_P422_interval"]]
        p447 = [float(value) for value in result["mse_P447_interval"]]
        self.assertLess(p447[1], p422[0])
        self.assertLess(p422[1], baseline[0])
        reduction = [float(value) for value in result["reduction_vs_P422_interval"]]
        self.assertGreater(reduction[0], 0.0034)

    def test_local_only_artifacts(self) -> None:
        for path in (
            batch.RESULTS_PATH, batch.SUMMARY_PATH, batch.P445_PATH,
            batch.P446_PATH, batch.P446_SIMPLEX_PATH, batch.P447_PATH,
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
