#!/usr/bin/env python3
"""Tests for the strict-alpha radical correction scratch probe."""
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_radical_correction_probe.py"
REPORT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_radical_correction_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_radical_correction_report.md"


class BridgeStrictAlphaRadicalCorrectionProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(REPORT.read_text(encoding="utf-8"))

    def test_result_kind_and_guardrails(self) -> None:
        self.assertEqual(self.payload["result_kind"], "SCRATCH_STRICT_ALPHA_RADICAL_CORRECTION_PROBE__NOT_A_THEOREM")
        self.assertIn("No identity K_legacy_ont == K_strict_gate is exported.", self.payload["hard_limits"])
        self.assertIn("No legacy physical-role transfer is licensed.", self.payload["hard_limits"])
        self.assertIn("No QW-2191 selector discharge and no ToE closure are claimed.", self.payload["hard_limits"])

    def test_radical_target_is_exact_for_eta(self) -> None:
        target = self.payload["exact_radical_target"]
        self.assertAlmostEqual(target["q_radical"], 1.010477711006932, places=15)
        self.assertLess(abs(target["polynomial_residual_243_q5_minus_256"]), 1e-12)
        self.assertAlmostEqual(target["eta_from_radical"], 1.8, places=15)
        self.assertLess(abs(target["eta_residual_vs_9_5"]), 1e-12)

    def test_rational_shadow_comparison(self) -> None:
        comparison = self.payload["comparison_to_rational_shadows"]
        self.assertFalse(comparison["diophantine_exact_rational_slot_possible"])
        self.assertLess(abs(comparison["natural_96_over_95"]["eta_residual_vs_9_5"]), 5e-5)
        self.assertLess(abs(comparison["best_bounded_diophantine_row"]["eta_signed_residual_vs_9_5"]), 1e-6)

    def test_candidate_supported_but_not_derivation(self) -> None:
        candidate = self.payload["candidate_interpretation"]
        self.assertTrue(candidate["supported_by_this_probe"])
        self.assertIn("not to any rational slot exclusion", candidate["content"])
        self.assertIn("No strict-side", candidate["why_this_is_not_enough"])
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
