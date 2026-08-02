#!/usr/bin/env python3
"""Regression tests for FIN Programs P458/P459/P464."""

import json
from fractions import Fraction
from pathlib import Path
import unittest

import fin_programs_458_459_464 as p

ROOT = Path(__file__).resolve().parent


class Programs458459464Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.results = json.loads((ROOT / "FIN_Programs_458_459_464_Results.json").read_text(encoding="utf-8"))

    def test_p458_unique_palindromic_maximizer(self) -> None:
        r = self.results["P458"]
        self.assertLess(r["maximum_second_derivative_upper"], -3.8)
        self.assertGreater(r["left_endpoint_derivative_interval"][0], 0)
        self.assertLess(r["right_endpoint_derivative_interval"][1], 0)
        self.assertLess(r["unique_maximizer_width"], 1.3e-13)
        self.assertIn("not claimed", r["boundary"])

    def test_p459_fiberwise_not_physical(self) -> None:
        r = self.results["P459"]
        self.assertTrue(r["fiberwise_canonical"])
        self.assertFalse(r["box_collapsed_to_single_physical_allocation"])
        self.assertGreater(float(r["minimum_certified_hessian_diagonal_lower"]), 90)
        self.assertGreater(float(r["minimum_adjacent_ranking_margin"]), 0.005)
        self.assertIn("synthetic", r["physical_boundary"])

    def test_p464_exact_improvement_and_slacks(self) -> None:
        r = self.results["P464"]
        self.assertTrue(r["strict_improvement_accepted"])
        self.assertTrue(r["candidate_rational_certificate_feasible"])
        self.assertLess(Fraction(r["accepted_global_gap"]), Fraction(23, 100_000_000))
        self.assertLess(Fraction(r["accepted_certified_upper"]), Fraction(r["previous_certified_upper"]))
        for row in r["candidate_certificate"]["slacks"]:
            self.assertGreater(float(Fraction(row["certified_exact_minimum_eigenvalue_lower"])), 0)

    def test_metadata_boundaries(self) -> None:
        m = self.results["metadata"]
        for key in ("external_physical_evidence", "selector_discharged", "dimensional_source_exported", "legacy_strict_bridge_complete", "legacy_role_transfer_started"):
            self.assertFalse(m[key])


if __name__ == "__main__":
    unittest.main()

