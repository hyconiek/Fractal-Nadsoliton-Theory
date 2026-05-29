#!/usr/bin/env python3
"""Tests for the strict-alpha fractal dimension scratch probe."""
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_fractal_dimension_probe.py"
REPORT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_fractal_dimension_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_fractal_dimension_report.md"


class BridgeStrictAlphaFractalDimensionProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(REPORT.read_text(encoding="utf-8"))

    def test_result_kind_and_guardrails(self) -> None:
        self.assertEqual(self.payload["result_kind"], "SCRATCH_STRICT_ALPHA_FRACTAL_DIMENSION_PROBE__NOT_A_THEOREM")
        self.assertIn("No identity K_legacy_ont == K_strict_gate is exported.", self.payload["hard_limits"])
        self.assertIn("No legacy physical-role transfer is licensed.", self.payload["hard_limits"])
        self.assertIn("No QW-2191 selector discharge and no ToE closure are claimed.", self.payload["hard_limits"])

    def test_alpha_z12_dimension_is_close_but_not_exact(self) -> None:
        candidate = self.payload["fractal_dimension_candidate"]
        self.assertAlmostEqual(candidate["eta_alpha_z12"], 1.7924812503605783, places=12)
        self.assertTrue(candidate["within_one_percent_of_9_5"])
        self.assertGreater(candidate["abs_residual_vs_9_5"], 1e-3)
        self.assertLess(candidate["relative_abs_residual_vs_9_5"], 0.01)

    def test_inverse_and_integer_scale_discriminators(self) -> None:
        inverse = self.payload["inverse_checks"]
        best = self.payload["integer_scale_discriminator"]["best_integer_scale_by_abs_residual"]
        self.assertEqual(best["scale"], 4)
        self.assertLess(inverse["relative_scale_residual"], 0.01)
        self.assertGreater(inverse["relative_effective_count_residual"], 0.005)
        self.assertLess(inverse["relative_effective_count_residual"], 0.02)

    def test_candidate_and_upstream_replays(self) -> None:
        interpretation = self.payload["candidate_interpretation"]
        replay = self.payload["upstream_replay"]
        self.assertTrue(interpretation["supported_by_this_probe"])
        self.assertTrue(replay["compression_audit_says_strict_terminal_candidate"])
        self.assertFalse(replay["compression_audit_says_not_final_formula"])
        self.assertTrue(replay["measure_transport_candidate_supported"])
        self.assertLess(replay["measure_transport_balance_residual"], 1e-12)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
