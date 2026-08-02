#!/usr/bin/env python3
"""Regression tests for local FIN Programs P454, P455, and P457."""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import unittest

import numpy as np

import fin_programs_454_455_457 as current


ROOT = Path(__file__).resolve().parent


class Programs454455457Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = json.loads(current.RESULTS_PATH.read_text(encoding="utf-8"))

    def test_p454_rational_dual_and_global_bracket(self) -> None:
        result = self.data["P454"]
        lower = Fraction(result["inherited_primal_lower"])
        upper = Fraction(result["certified_dual_upper"])
        gap = Fraction(result["certified_global_gap"])
        self.assertGreater(upper, lower)
        self.assertEqual(upper - lower, gap)
        self.assertLess(float(gap), 6.7e-6)
        self.assertTrue(result["rational_dual_certificate"]["feasible"])
        self.assertEqual(len(result["rational_dual_certificate"]["slacks"]), 8)
        for slack in result["rational_dual_certificate"]["slacks"]:
            self.assertGreater(
                float(Fraction(slack["certified_exact_minimum_eigenvalue_lower"])), 0
            )
        payload = np.load(current.P454_WITNESS_PATH)
        self.assertEqual(payload["X3"].shape, (8, 8))
        self.assertEqual(payload["X2"].shape, (4, 4))
        self.assertEqual(payload["X1"].shape, (2, 2))
        self.assertAlmostEqual(float(payload["lambda"][0]), float(upper), places=12)

    def test_p455_exact_symmetry_residual_and_numerical_falsification(self) -> None:
        result = self.data["P455"]
        self.assertEqual(result["real_c3_affine_dimension"], 14)
        self.assertEqual(result["complement_constraint_rank"], 9)
        self.assertEqual(result["real_complement_fixed_dimension"], 5)
        self.assertEqual(result["O167_face_dimension"], 3)
        self.assertEqual(result["symmetry_allowed_residual_dimensions"], 2)
        self.assertTrue(result["face_contained_in_fixed_space"])
        self.assertEqual(result["fixed_space_starts"], 8)
        self.assertLess(abs(result["maximum_fixed_minus_face_gain"]), 1e-12)
        self.assertLess(max(result["fixed_space_hessian_eigenvalues"]), -0.3)
        self.assertIn("cannot by themselves", result["obstruction"])

    def test_p457_refined_global_cover(self) -> None:
        result = self.data["P457"]
        self.assertGreater(result["global_upper"], result["global_lower"])
        self.assertLess(result["global_gap"], 1e-6)
        self.assertEqual(result["evaluated_boxes"], 5008)
        self.assertEqual(result["terminal_boxes"], 5009)
        self.assertGreater(result["improvement_factor_over_P452_gap"], 990)
        self.assertLess(result["maximizer_hull"][0], result["incumbent_a"])
        self.assertGreater(result["maximizer_hull"][1], result["incumbent_a"])

    def test_artifacts_and_epistemic_boundaries(self) -> None:
        for path in (
            current.RESULTS_PATH,
            current.SUMMARY_PATH,
            current.P454_PATH,
            current.P454_WITNESS_PATH,
            current.P455_PATH,
            current.P457_PATH,
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
