#!/usr/bin/env python3
"""Regression tests for FIN local Programs P471/P472/P473."""

import json
from fractions import Fraction
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parent


class Programs471472473Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.results = json.loads(
            (ROOT / "FIN_Programs_471_472_473_Results.json").read_text(encoding="utf-8")
        )

    def test_p471_exact_polynomial_reduction(self) -> None:
        result = self.results["P471"]
        self.assertEqual(result["upper_triangle_equations"], 36)
        self.assertEqual(result["exact_residual_orbits"], 13)
        self.assertEqual(result["polynomial_degree_maximum"], 3)
        self.assertLess(result["floating_residual_infinity_norm"], 1e-15)
        self.assertGreater(result["jacobian_singular_values"][-1], 0.003)
        self.assertIn("do not prove", result["boundary"])

    def test_p472_exact_krawczyk_inclusion(self) -> None:
        result = self.results["P472"]
        self.assertTrue(result["exact_rational_interval_arithmetic"])
        self.assertTrue(result["strict_inclusion"])
        self.assertEqual(Fraction(result["selected_radius"]), Fraction(3, 10**14))
        self.assertGreater(Fraction(result["minimum_inclusion_margin"]), 0)
        self.assertNotIn("1/100000000000000", result["admitted_radii"])

    def test_p473_exact_attainment_but_not_global_uniqueness(self) -> None:
        result = self.results["P473"]
        self.assertTrue(result["exact_root_proved"])
        self.assertTrue(result["root_unique_in_box"])
        self.assertTrue(result["exact_O167_attainment_proved"])
        self.assertFalse(result["full_cone_optimizer_uniqueness_proved"])
        self.assertEqual(
            Fraction(result["exact_global_value_interval_width"]),
            Fraction(3, 5 * 10**13),
        )
        self.assertGreater(Fraction(result["normalizer_box_positive_lower"]), 0)
        self.assertGreater(Fraction(result["x3_box_positive_lower"]), 0)
        self.assertIn("does not by itself exclude", result["boundary"])

    def test_metadata_boundaries(self) -> None:
        metadata = self.results["metadata"]
        for key in (
            "external_physical_evidence",
            "selector_discharged",
            "dimensional_source_exported",
            "legacy_strict_bridge_complete",
            "legacy_role_transfer_started",
        ):
            self.assertFalse(metadata[key])


if __name__ == "__main__":
    unittest.main()
