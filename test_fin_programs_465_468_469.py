#!/usr/bin/env python3
"""Regression tests for FIN local Programs P465/P468/P469."""

import json
from fractions import Fraction
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parent


class Programs465468469Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.results = json.loads(
            (ROOT / "FIN_Programs_465_468_469_Results.json").read_text(encoding="utf-8")
        )

    def test_p465_strict_full_simplex_uniqueness(self) -> None:
        result = self.results["P465"]
        self.assertIn("complete four-sector simplex", result["status"])
        self.assertGreater(
            float(Fraction(result["strict_condition_interval_kappa_minus_one"][0])),
            1.13,
        )
        root = result["full_simplex_unique_maximizer_interval"]
        self.assertLess(root[1] - root[0], 1.3e-13)
        self.assertIn("not an unrestricted", result["proof_boundary"])

    def test_p468_exact_rational_gap(self) -> None:
        result = self.results["P468"]
        self.assertTrue(result["strict_improvement_accepted"])
        self.assertLess(
            Fraction(result["accepted_global_gap"]), Fraction(1, 100_000_000)
        )
        self.assertLess(
            Fraction(result["accepted_certified_upper"]),
            Fraction(result["previous_certified_upper"]),
        )
        self.assertGreater(result["gap_reduction_factor_over_P464"], 500)
        for row in result["accepted_certificate"]["slacks"]:
            self.assertGreater(
                float(Fraction(row["certified_exact_minimum_eigenvalue_lower"])), 0
            )
        self.assertIn("does not prove exact", result["boundary"])

    def test_p469_is_strong_evidence_not_exact_proof(self) -> None:
        result = self.results["P469"]
        self.assertIn("[Strong evidence]", result["status"])
        self.assertIn("[Open]", result["status"])
        self.assertLess(result["full_residual_infinity_norm"], 2e-15)
        self.assertGreater(result["normalizer_minimum_eigenvalue"], 0.05)
        self.assertGreater(result["transformed_delta_minimum_absolute_eigenvalue"], 0.008)
        self.assertIn("Krawczyk", result["boundary"])

    def test_metadata_boundaries(self) -> None:
        metadata = self.results["metadata"]
        for key in (
            "external_physical_evidence",
            "selector_discharged",
            "dimensional_source_exported",
            "legacy_strict_bridge_complete",
            "legacy_role_transfer_started",
            "exact_O167_attainment_proven",
        ):
            self.assertFalse(metadata[key])


if __name__ == "__main__":
    unittest.main()
