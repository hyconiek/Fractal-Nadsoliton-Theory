#!/usr/bin/env python3
"""Tests for the strict-alpha finite-size correction scratch probe."""
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_finite_size_correction_probe.py"
REPORT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_finite_size_correction_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_finite_size_correction_report.md"


class BridgeStrictAlphaFiniteSizeCorrectionProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(REPORT.read_text(encoding="utf-8"))

    def test_result_kind_and_guardrails(self) -> None:
        self.assertEqual(self.payload["result_kind"], "SCRATCH_STRICT_ALPHA_FINITE_SIZE_CORRECTION_PROBE__NOT_A_THEOREM")
        self.assertIn("No identity K_legacy_ont == K_strict_gate is exported.", self.payload["hard_limits"])
        self.assertIn("No legacy physical-role transfer is licensed.", self.payload["hard_limits"])
        self.assertIn("No QW-2191 selector discharge and no ToE closure are claimed.", self.payload["hard_limits"])

    def test_exact_correction_target(self) -> None:
        exact = self.payload["exact_correction_required"]
        self.assertAlmostEqual(exact["q_target"], 1.010477711006932, places=15)
        self.assertAlmostEqual(exact["effective_count_target_12_times_q"], 12.125732532083182, places=12)
        self.assertGreater(exact["relative_count_correction"], 0.01)
        self.assertLess(exact["relative_count_correction"], 0.011)

    def test_bounded_rational_search(self) -> None:
        search = self.payload["bounded_rational_search"]
        best = search["best_candidate"]
        self.assertEqual(best["reduced_label"], "96/95")
        self.assertLess(best["eta_abs_residual_vs_9_5"], 5e-5)
        self.assertGreater(search["baseline_to_best_eta_residual_improvement_factor"], 100.0)
        self.assertEqual(search["nearest_fraction_limit_denominator"]["label"], "96/95")

    def test_candidate_is_supported_but_not_theorem(self) -> None:
        candidate = self.payload["candidate_interpretation"]
        self.assertTrue(candidate["supported_by_this_probe"])
        self.assertIn("no strict-side", candidate["why_this_is_not_enough"].lower())
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
