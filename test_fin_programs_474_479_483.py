#!/usr/bin/env python3
"""Regression tests for FIN local Programs P474, P479, and P483."""

from fractions import Fraction
import json
from pathlib import Path
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "FIN_Programs_474_479_483_Results.json"


class Programs474479483Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.results = json.loads(RESULTS.read_text(encoding="utf-8"))

    def test_p474_is_not_overpromoted(self) -> None:
        result = self.results["P474"]
        self.assertIn("[Strong evidence]", result["status"])
        self.assertIn("[Open]", result["status"])
        self.assertEqual(result["full_contact_nullity"], 1)
        self.assertEqual(result["imaginary_causal_riccati_nullity"], 1)
        self.assertGreater(result["smallest_nonzero_contact_singular_value"], 1e-3)
        self.assertLess(result["smallest_contact_singular_value"], 1e-12)
        self.assertGreater(result["minimum_scanned_normalizer_eigenvalue"], 0)
        self.assertLess(result["maximum_scanned_objective_deviation"], 2e-14)

    def test_p474_witness_residuals(self) -> None:
        witness = np.load(ROOT / "FIN_Programs_474_479_483_P474_Flat_Direction.npz")
        x3 = witness["X3"]
        delta = witness["Delta"]
        q_direction = witness["Q_direction"]
        k_matrix = np.real(delta / 1j)
        residual = x3 @ q_direction @ x3 + k_matrix @ q_direction @ k_matrix / 4
        self.assertLess(np.linalg.norm(residual, np.inf), 2e-13)
        self.assertLess(np.linalg.norm(q_direction + q_direction.T), 2e-12)

    def test_p479_records_local_checker_boundary(self) -> None:
        result = self.results["P479"]
        self.assertFalse(result["lean_checked"])
        self.assertTrue(result["local_toolchain_missing"])
        self.assertIn("not described as machine-checked", result["boundary"])
        source = ROOT / result["formal_source"]
        self.assertTrue(source.exists())
        text = source.read_text(encoding="utf-8")
        self.assertIn("exact_attainment_blocks_strict_improvement", text)

    def test_p483_uniform_exact_tube(self) -> None:
        result = self.results["P483"]
        self.assertTrue(result["uniform_root_tube_proved"])
        self.assertGreaterEqual(result["admitted_trials"], 1)
        selected = result["selected_certificate"]
        self.assertTrue(selected["included"])
        self.assertGreater(Fraction(selected["minimum_inclusion_margin"]), 0)
        self.assertLess(Fraction(selected["maximum_contraction_row_sum"]), 1)
        self.assertEqual(Fraction(selected["q_radius"]), Fraction(3, 10**9))
        self.assertEqual(Fraction(selected["theta_radius"]), Fraction(3, 10**9))
        self.assertGreater(Fraction(result["normalizer_tube_positive_lower"]), 0)
        self.assertGreater(Fraction(result["X3_tube_positive_lower"]), 0)

    def test_metadata_boundaries(self) -> None:
        metadata = self.results["metadata"]
        self.assertFalse(metadata["network_used"])
        self.assertFalse(metadata["laboratory_data_used"])
        self.assertFalse(metadata["external_audit_used"])
        self.assertIn("QW-2191 remains open", metadata["selector_boundary"])
        self.assertIn("No legacy/strict substitution", metadata["kernel_boundary"])


if __name__ == "__main__":
    unittest.main()
