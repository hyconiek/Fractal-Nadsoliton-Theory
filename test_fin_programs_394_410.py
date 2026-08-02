#!/usr/bin/env python3
"""Regression tests for FIN Programs P394--P410 (Release 10.34)."""

from __future__ import annotations

import csv
from fractions import Fraction
import json
from pathlib import Path
import unittest

import fin_p403_moment_validator as validator


ROOT = Path(__file__).resolve().parent
RESULTS = json.loads((ROOT / "FIN_Programs_394_410_Results.json").read_text(encoding="utf-8"))


class Programs394To410Tests(unittest.TestCase):
    def test_all_programs_reported(self) -> None:
        self.assertEqual(
            [key for key in RESULTS if key.startswith("P")],
            [f"P{index}" for index in range(394, 411)],
        )

    def test_p394_exact_root_brackets(self) -> None:
        with (ROOT / "FIN_Programs_394_410_Trust_Boundary.csv").open(
            encoding="utf-8", newline=""
        ) as handle:
            rows = list(csv.DictReader(handle))[:11]
        self.assertEqual(len(rows), 11)
        for row in rows:
            d = int(row["distance"])
            lo = Fraction(row["lower_exact"])
            hi = Fraction(row["upper_exact"])
            self.assertLessEqual(lo**5, d)
            self.assertLess(d, hi**5)
            self.assertEqual(hi - lo, Fraction(1, 10**18))
        self.assertEqual(RESULTS["P394"]["lean_returncode"], 0)
        self.assertEqual(RESULTS["P394"]["adaptive_nodes_visited"], 147)
        self.assertEqual(RESULTS["P394"]["maximum_uniform_cells_covered"], 2**14)

    def test_p395_exact_root_count_and_slack(self) -> None:
        result = RESULTS["P395"]["frozen_P380_exact_Sturm_audit"]
        self.assertEqual(result["exact_derivative_roots_in_unit_interval"], 9)
        self.assertEqual(result["isolated_root_boxes"], 9)
        self.assertLessEqual(result["maximum_root_box_width"], 2**-50)
        self.assertEqual(result["positive_primal_slack_atoms"], 7)
        self.assertGreater(result["minimum_primal_slack_lower"], 0)
        self.assertLess(RESULTS["P395"]["maximum_80_digit_residual"], 1e-25)

    def test_p397_exact_integer_counterexample(self) -> None:
        defect = 30599**5 - 512 * 10201**5
        self.assertEqual(defect, RESULTS["P397"]["exact_fifth_power_integer_defect"])
        self.assertNotEqual(defect, 0)

    def test_p398_bounded_recovery(self) -> None:
        result = RESULTS["P398"]
        self.assertEqual(result["multi_starts"], 4)
        self.assertEqual(result["optimizer_successes"], 4)
        self.assertEqual(result["distant_alias_candidates"], 0)

    def test_p400_symmetric_codebook_improves_baselines(self) -> None:
        result = RESULTS["P400"]
        self.assertEqual(result["grid_rows"], 64)
        self.assertGreater(result["maximum_gain_over_product_or_ghz"], 0.05)
        self.assertGreater(result["maximum_remaining_adaptive_gap"], 0.4)

    def test_p401_extremal_two_mode_survives_full_mode_search(self) -> None:
        result = RESULTS["P401"]
        self.assertLess(result["fourier_off_diagonal_defect"], 1e-12)
        self.assertLess(result["maximum_gain_over_extremal_two_mode"], 1e-12)

    def test_p402_is_finite_grid_maximin(self) -> None:
        result = RESULTS["P402"]
        self.assertEqual(result["design_count"], 12)
        self.assertEqual(result["grid_points_per_design"], 1501)
        self.assertGreater(result["best_worst_case_value"], 0.99)

    def test_p403_validator_and_boundary(self) -> None:
        result = validator.validate(
            ROOT / "FIN_P403_Jordan_Sampling_Protocol.json",
            ROOT / "FIN_P403_Jordan_Sampling_Synthetic_Events.csv",
        )
        self.assertTrue(result["validation_pass"])
        self.assertTrue(result["synthetic_reference"])
        self.assertFalse(result["physical_evidence_admitted"])
        self.assertEqual(result["event_count"], 50000)
        self.assertLess(result["maximum_moment_z_score"], 6)

    def test_p405_selector_not_discharged(self) -> None:
        self.assertFalse(RESULTS["P405"]["selector_discharge"])

    def test_p406_conditional_scale_only(self) -> None:
        result = RESULTS["P406"]
        self.assertAlmostEqual(result["rho_midpoint"], 0.4296970877901551)
        self.assertLess(result["enclosure_width"], 1e-50)
        self.assertIn("added normalization", result["boundary"].lower())
        self.assertIn("does not create", result["boundary"].lower())

    def test_external_gates_remain_empty(self) -> None:
        for program in range(407, 411):
            result = RESULTS[f"P{program}"]
            self.assertFalse(result["admitted"])

    def test_pdf_exists_and_nonempty(self) -> None:
        pdf = ROOT / "FIN_Programs_394_410_Formal_Operational_Report_EN.pdf"
        self.assertTrue(pdf.exists())
        self.assertGreater(pdf.stat().st_size, 300_000)


if __name__ == "__main__":
    unittest.main()
