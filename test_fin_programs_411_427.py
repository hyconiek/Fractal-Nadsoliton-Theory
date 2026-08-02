#!/usr/bin/env python3
"""Regression tests for FIN Programs P411--P427."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import subprocess
import unittest


ROOT = Path(__file__).resolve().parent
RESULTS = json.loads((ROOT / "FIN_Programs_411_427_Results.json").read_text(encoding="utf-8"))
LEAN = ROOT / ".elan/toolchains/leanprover--lean4---v4.28.0/bin/lean"


def table(name: str) -> list[dict[str, str]]:
    with (ROOT / name).open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


class Programs411To427Tests(unittest.TestCase):
    def test_all_programs_reported(self) -> None:
        self.assertEqual([key for key in RESULTS if key.startswith("P")], [f"P{x}" for x in range(411, 428)])

    def test_p411_lean_and_widths(self) -> None:
        self.assertEqual(RESULTS["P411"]["lean_returncode"], 0)
        self.assertLess(RESULTS["P411"]["maximum_rational_width"], 1e-30)
        completed = subprocess.run([str(LEAN), "FIN_Programs_411_427_Taylor_Cosine.lean"], cwd=ROOT, capture_output=True, text=True, timeout=120)
        self.assertEqual(completed.returncode, 0, completed.stderr)

    def test_p412_contact_candidate_is_not_overpromoted(self) -> None:
        self.assertLess(RESULTS["P412"]["maximum_80_digit_residual"], 1e-25)
        self.assertIn("[Open]", RESULTS["P412"]["status"])

    def test_p413_names_transcendence_gap(self) -> None:
        self.assertIn("Lindemann-Weierstrass", RESULTS["P413"]["formal_gap"])

    def test_p414_exact_obstruction(self) -> None:
        self.assertNotEqual(RESULTS["P414"]["exact_integer_defect"], 0)
        self.assertIn("[Refuted]", RESULTS["P414"]["status"])

    def test_p415_sampled_full_rank(self) -> None:
        self.assertEqual(RESULTS["P415"]["full_rank_samples"], RESULTS["P415"]["samples"])
        self.assertIn("[Open]", RESULTS["P415"]["status"])

    def test_external_gates_are_closed(self) -> None:
        for program in (416, 420, 421, 426, 427):
            self.assertFalse(RESULTS[f"P{program}"]["admitted"])

    def test_p417_bounds_are_ordered(self) -> None:
        for row in table("FIN_Programs_411_427_Noisy_Comb_Gap.csv"):
            self.assertLessEqual(float(row["optimized_parallel_lower"]), float(row["adaptive_hybrid_upper"]) + 2e-8)

    def test_p418_probabilities_and_gains(self) -> None:
        for row in table("FIN_Programs_411_427_Erasure_Robust_Code.csv"):
            probabilities = json.loads(row["sector_probabilities"])
            self.assertTrue(math.isclose(sum(probabilities), 1.0, abs_tol=1e-9))
            self.assertGreaterEqual(float(row["gain_over_product_or_ghz"]), -2e-8)

    def test_p419_fourier_diagonal(self) -> None:
        self.assertLess(RESULTS["P419"]["fourier_off_diagonal_defect"], 1e-10)

    def test_p422_variance_is_reduced(self) -> None:
        self.assertGreater(RESULTS["P422"]["relative_variance_reduction"], 0)
        self.assertLess(RESULTS["P422"]["optimal_total_variance"], RESULTS["P422"]["baseline_total_variance"])

    def test_p423_does_not_claim_selector(self) -> None:
        self.assertFalse(RESULTS["P423"]["selector_discharge"])
        self.assertIn("QW-2191", RESULTS["P423"]["boundary"])

    def test_p424_is_odd_but_not_selector(self) -> None:
        self.assertTrue(RESULTS["P424"]["nonzero"])
        self.assertLess(abs(RESULTS["P424"]["oddness_defect"]), 1e-12)
        self.assertFalse(RESULTS["P424"]["selector_discharge"])

    def test_p425_sections_are_distinct(self) -> None:
        self.assertTrue(RESULTS["P425"]["pairwise_distinct_baseline_elasticities"])

    def test_summary_and_figures(self) -> None:
        self.assertEqual(len(table("FIN_Programs_411_427_Summary.csv")), 17)
        expected = {"p411_p412_formal_contact.png", "p414_p415_damping_photonic.png", "p417_p418_noisy_erasure.png", "p422_p425_estimator_scale.png"}
        actual = {path.name for path in (ROOT / "FIN_Programs_411_427_Figures").glob("*.png")}
        self.assertEqual(actual, expected)
        for name in expected:
            self.assertGreater((ROOT / "FIN_Programs_411_427_Figures" / name).stat().st_size, 10_000)

    def test_json_is_strict(self) -> None:
        text = (ROOT / "FIN_Programs_411_427_Results.json").read_text(encoding="utf-8")
        self.assertNotIn("NaN", text)
        self.assertNotIn("Infinity", text)


if __name__ == "__main__":
    unittest.main()
