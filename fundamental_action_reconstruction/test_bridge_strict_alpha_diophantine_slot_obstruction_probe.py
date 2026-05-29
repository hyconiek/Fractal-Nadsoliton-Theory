#!/usr/bin/env python3
"""Tests for the strict-alpha Diophantine slot obstruction scratch probe."""
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_diophantine_slot_obstruction_probe.py"
REPORT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_diophantine_slot_obstruction_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_diophantine_slot_obstruction_report.md"


class BridgeStrictAlphaDiophantineSlotObstructionProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(REPORT.read_text(encoding="utf-8"))

    def test_result_kind_and_guardrails(self) -> None:
        self.assertEqual(
            self.payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_DIOPHANTINE_SLOT_OBSTRUCTION_PROBE__NOT_A_THEOREM",
        )
        self.assertIn("No identity K_legacy_ont == K_strict_gate is exported.", self.payload["hard_limits"])
        self.assertIn("No legacy physical-role transfer is licensed.", self.payload["hard_limits"])
        self.assertIn("No QW-2191 selector discharge and no ToE closure are claimed.", self.payload["hard_limits"])

    def test_exact_rational_slot_exclusion_is_obstructed(self) -> None:
        obstruction = self.payload["exact_obstruction"]
        self.assertFalse(obstruction["exact_rational_slot_exclusion_possible"])
        self.assertIn("irrational", obstruction["irrationality_reason"])
        self.assertAlmostEqual(obstruction["q_target"], 1.010477711006932, places=15)
        self.assertAlmostEqual(obstruction["exact_B_over_k_if_formally_solving_B_over_B_minus_k"], 96.44069304244016, places=12)

    def test_natural_shadow_and_bounded_search_are_separated(self) -> None:
        natural = self.payload["natural_96_shadow_replay"]["natural_one_slot_row"]
        search = self.payload["bounded_integer_search"]
        best = search["best_row"]
        self.assertEqual(natural["correction_label"], "96/95")
        self.assertLess(natural["eta_abs_residual_vs_9_5"], 5e-5)
        self.assertFalse(search["best_is_natural_96_one_slot"])
        self.assertEqual(best["slot_count_B"], 1543)
        self.assertEqual(best["excluded_slots_k"], 16)
        self.assertLess(best["eta_abs_residual_vs_9_5"], 1e-6)

    def test_candidate_supported_but_downgraded_to_approximants(self) -> None:
        candidate = self.payload["candidate_interpretation"]
        self.assertTrue(candidate["supported_by_this_probe"])
        self.assertIn("obstructed", candidate["content"])
        self.assertIn("approximant", self.payload["bounded_integer_search"]["natural_96_rank_message"])
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
