#!/usr/bin/env python3
"""Tests for the strict-alpha balanced branch ledger scratch probe."""
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_balanced_branch_ledger_probe.py"
REPORT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_balanced_branch_ledger_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_balanced_branch_ledger_report.md"


class BridgeStrictAlphaBalancedBranchLedgerProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(REPORT.read_text(encoding="utf-8"))

    def test_result_kind_and_guardrails(self) -> None:
        self.assertEqual(self.payload["result_kind"], "SCRATCH_STRICT_ALPHA_BALANCED_BRANCH_LEDGER_PROBE__NOT_A_THEOREM")
        self.assertIn("No identity K_legacy_ont == K_strict_gate is exported.", self.payload["hard_limits"])
        self.assertIn("No legacy physical-role transfer is licensed.", self.payload["hard_limits"])
        self.assertIn("No theorem derives the balance/equipartition selector as a strict-core selector source.", self.payload["hard_limits"])
        self.assertIn("No QW-2191 selector discharge and no ToE closure are claimed.", self.payload["hard_limits"])

    def test_weakened_assumption_replaces_one_or_two_bit_axiom(self) -> None:
        weakened = self.payload["weakened_assumptions_compared_to_previous_probe"]
        self.assertIn("{1,2}", weakened["previous_extra_assumption"])
        self.assertIn("positive integer", weakened["replacement_assumption"])
        self.assertIn("balance selector", weakened["replacement_assumption"])
        self.assertIn("forced", weakened["forced_consequence"])

    def test_positive_integer_scan_selects_unique_balanced_ledger(self) -> None:
        scan = self.payload["positive_integer_ledger_scan"]
        self.assertEqual(scan["labelled_positive_ledgers_count"], 35)
        self.assertEqual(scan["canonical_positive_ledgers_count"], 3)
        self.assertEqual(scan["variance_minimizers"], [[2, 2, 2, 1, 1]])
        self.assertTrue(scan["variance_minimizer_unique"])
        self.assertEqual(scan["max_gap_minimizers"], [[2, 2, 2, 1, 1]])
        self.assertTrue(scan["max_gap_minimizer_unique"])
        self.assertEqual(scan["min_variance_score_fraction"], "6/25")

    def test_selected_balanced_ledger_exact_replay(self) -> None:
        selected = self.payload["selected_balanced_ledger"]
        self.assertEqual(selected["ledger"], [2, 2, 2, 1, 1])
        self.assertEqual(selected["branch_ratio_labels"], ["4/3", "4/3", "4/3", "2/3", "2/3"])
        self.assertEqual(selected["product_fraction"], "256/243")
        self.assertLess(abs(selected["geometric_mean_residual_vs_q_radical"]), 1e-15)
        self.assertLess(abs(selected["eta_residual_vs_9_5"]), 1e-12)

    def test_candidate_supported_but_selector_not_derived(self) -> None:
        candidate = self.payload["candidate_interpretation"]
        self.assertTrue(candidate["supported_by_this_probe"])
        self.assertIn("weaker positive-integer plus balance-selector premise", candidate["content"])
        self.assertIn("not derived", candidate["why_this_is_not_enough"])
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
