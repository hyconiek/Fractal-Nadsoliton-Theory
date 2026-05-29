#!/usr/bin/env python3
"""Tests for the strict-alpha binary/ternary ledger scratch probe."""
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_binary_ternary_ledger_probe.py"
REPORT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_binary_ternary_ledger_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_binary_ternary_ledger_report.md"


class BridgeStrictAlphaBinaryTernaryLedgerProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(REPORT.read_text(encoding="utf-8"))

    def test_result_kind_and_guardrails(self) -> None:
        self.assertEqual(self.payload["result_kind"], "SCRATCH_STRICT_ALPHA_BINARY_TERNARY_LEDGER_PROBE__NOT_A_THEOREM")
        self.assertIn("No identity K_legacy_ont == K_strict_gate is exported.", self.payload["hard_limits"])
        self.assertIn("No legacy physical-role transfer is licensed.", self.payload["hard_limits"])
        self.assertIn("No QW-2191 selector discharge and no ToE closure are claimed.", self.payload["hard_limits"])

    def test_exact_target_prime_ledger_replay(self) -> None:
        replay = self.payload["exact_target_replay"]
        self.assertAlmostEqual(replay["q_radical"], 1.010477711006932, places=15)
        self.assertAlmostEqual(replay["q_radical_power_5"], 256 / 243, places=15)
        self.assertEqual(replay["target_binary_exponent_total"], 8)
        self.assertEqual(replay["target_ternary_exponent_total"], -5)

    def test_unique_one_or_two_bit_ledger(self) -> None:
        solution = self.payload["ledger_solution_under_candidate_rule"]
        self.assertTrue(solution["unique_under_rule"])
        self.assertEqual(solution["solutions_as_canonical_multisets"], [[2, 2, 2, 1, 1]])
        self.assertEqual(solution["selected_count_by_exponent"], {"2": 3, "1": 2})
        self.assertEqual([row["branch_ratio_label"] for row in solution["branch_rows"]], ["4/3", "4/3", "4/3", "2/3", "2/3"])
        self.assertAlmostEqual(solution["product_branch_ratios"], 256 / 243, places=15)
        self.assertLess(abs(solution["geometric_mean_residual_vs_q_radical"]), 1e-15)
        self.assertLess(abs(solution["eta_residual_vs_9_5"]), 1e-12)

    def test_loose_ledger_warning_and_candidate_status(self) -> None:
        loose = self.payload["discriminator_against_looser_ledgers"]
        self.assertGreater(loose["number_of_canonical_ledgers_sum_8_length_5_with_entries_0_to_4"], 1)
        self.assertIn("branch law", loose["why_candidate_rule_matters"])
        candidate = self.payload["candidate_interpretation"]
        self.assertTrue(candidate["supported_by_this_probe"])
        self.assertIn("unique one-or-two-bit ledger", candidate["content"])
        self.assertIn("not derived", candidate["why_this_is_not_enough"])
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
