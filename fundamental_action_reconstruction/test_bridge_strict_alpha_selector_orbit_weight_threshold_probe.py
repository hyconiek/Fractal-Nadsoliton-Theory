#!/usr/bin/env python3
"""Tests for the strict-alpha selector orbit-weight threshold scratch probe."""
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_selector_orbit_weight_threshold_probe.py"
REPORT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_selector_orbit_weight_threshold_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_selector_orbit_weight_threshold_report.md"


class BridgeStrictAlphaSelectorOrbitWeightThresholdProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(REPORT.read_text(encoding="utf-8"))

    def test_result_kind_and_guardrails(self) -> None:
        self.assertEqual(self.payload["result_kind"], "SCRATCH_STRICT_ALPHA_SELECTOR_ORBIT_WEIGHT_THRESHOLD_PROBE__NOT_A_THEOREM")
        self.assertIn("No identity K_legacy_ont == K_strict_gate is exported.", self.payload["hard_limits"])
        self.assertIn("No legacy physical-role transfer is licensed.", self.payload["hard_limits"])
        self.assertIn("No theorem derives an intrinsic orbit-weight exponent gamma below log(3/2)/log(2).", self.payload["hard_limits"])
        self.assertIn("No QW-2191 selector discharge and no ToE closure are claimed.", self.payload["hard_limits"])

    def test_candidate_rows_share_eta_product(self) -> None:
        for row in self.payload["candidate_rows"]:
            self.assertEqual(row["product_fraction"], "256/243")
            self.assertLess(abs(row["eta_residual_vs_9_5"]), 1e-12)

    def test_critical_threshold_exact_inputs(self) -> None:
        threshold = self.payload["critical_threshold"]
        self.assertEqual(threshold["balanced_ledger_B"], [2, 2, 2, 1, 1])
        self.assertEqual(threshold["orbit_favoured_ledger_C"], [3, 2, 1, 1, 1])
        self.assertEqual(threshold["W_B_over_W_C"], "3/2")
        self.assertEqual(threshold["O_C_over_O_B"], "2/1")
        self.assertEqual(threshold["gamma_critical_closed_form"], "log(3/2)/log(2)")
        self.assertAlmostEqual(threshold["gamma_critical_numeric"], 0.5849625007211562, places=15)
        self.assertLess(abs(threshold["score_gap_B_minus_C_at_gamma_critical"]), 1e-12)

    def test_winner_scan_flips_at_threshold(self) -> None:
        scan = self.payload["winner_scan"]
        threshold = self.payload["critical_threshold"]
        self.assertEqual(scan["winners_by_gamma"]["0.000000"], [2, 2, 2, 1, 1])
        self.assertEqual(scan["winners_by_gamma"]["1.000000"], [3, 2, 1, 1, 1])
        self.assertEqual(threshold["winner_just_below_threshold"], [2, 2, 2, 1, 1])
        self.assertEqual(threshold["winner_just_above_threshold"], [3, 2, 1, 1, 1])
        self.assertTrue(scan["tie_at_gamma_c_between_B_and_C"])

    def test_candidate_supported_but_threshold_conditional(self) -> None:
        candidate = self.payload["candidate_interpretation"]
        self.assertTrue(candidate["supported_by_this_probe"])
        self.assertIn("below the exact orbit-weight threshold", candidate["content"])
        self.assertIn("exact bifurcation threshold", candidate["why_this_is_more_proof_like"])
        self.assertIn("No strict theorem", candidate["why_this_is_not_enough"])
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
