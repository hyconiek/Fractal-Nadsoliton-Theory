#!/usr/bin/env python3
"""Tests for the strict-alpha orbit/stabilizer correction scratch probe."""
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_orbit_stabilizer_correction_probe.py"
REPORT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_orbit_stabilizer_correction_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_orbit_stabilizer_correction_report.md"


class BridgeStrictAlphaOrbitStabilizerCorrectionProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(REPORT.read_text(encoding="utf-8"))

    def test_result_kind_and_guardrails(self) -> None:
        self.assertEqual(
            self.payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_ORBIT_STABILIZER_CORRECTION_PROBE__NOT_A_THEOREM",
        )
        self.assertIn("No identity K_legacy_ont == K_strict_gate is exported.", self.payload["hard_limits"])
        self.assertIn("No legacy physical-role transfer is licensed.", self.payload["hard_limits"])
        self.assertIn("No QW-2191 selector discharge and no ToE closure are claimed.", self.payload["hard_limits"])

    def test_natural_slot_count_and_one_slot_exclusion(self) -> None:
        counts = self.payload["strict_counts_used"]
        exact = self.payload["exact_vs_natural_one_slot_exclusion"]
        self.assertEqual(counts["natural_slot_count_B"], 96)
        self.assertAlmostEqual(exact["natural_q_B_over_B_minus_1"], 96 / 95, places=15)
        self.assertLess(exact["natural_eta_abs_residual_vs_9_5"], 5e-5)
        self.assertLess(abs(exact["natural_q_signed_residual_vs_exact"]), 5e-5)

    def test_slot_count_discriminator(self) -> None:
        discriminator = self.payload["slot_count_discriminator"]
        best = discriminator["best_natural_row"]
        self.assertTrue(discriminator["natural_96_is_best_in_candidate_library"])
        self.assertEqual(best["slot_count_B"], 96)
        self.assertEqual(best["correction_label"], "96/95")
        self.assertLess(best["eta_abs_residual_vs_9_5"], 5e-5)

    def test_candidate_supported_but_not_theorem(self) -> None:
        candidate = self.payload["candidate_interpretation"]
        self.assertTrue(candidate["supported_by_this_probe"])
        self.assertIn("No theorem", candidate["why_this_is_not_enough"])
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
