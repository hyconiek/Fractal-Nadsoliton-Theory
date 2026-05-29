#!/usr/bin/env python3
"""Tests for the strict-alpha minimal branch-count obstruction probe."""
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_minimal_branch_count_obstruction_probe.py"
REPORT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_minimal_branch_count_obstruction_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_minimal_branch_count_obstruction_report.md"


class BridgeStrictAlphaMinimalBranchCountObstructionProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(REPORT.read_text(encoding="utf-8"))

    def test_result_kind_and_guardrails(self) -> None:
        self.assertEqual(self.payload["result_kind"], "SCRATCH_STRICT_ALPHA_MINIMAL_BRANCH_COUNT_OBSTRUCTION_PROBE__NOT_A_THEOREM")
        self.assertIn("No identity K_legacy_ont == K_strict_gate is exported.", self.payload["hard_limits"])
        self.assertIn("No legacy physical-role transfer is licensed.", self.payload["hard_limits"])
        self.assertIn("No theorem derives physical five-branch channels; five is only the irreducible arithmetic count inside the stated model.", self.payload["hard_limits"])
        self.assertIn("No QW-2191 selector discharge and no ToE closure are claimed.", self.payload["hard_limits"])

    def test_exact_arithmetic_statement(self) -> None:
        statement = self.payload["exact_arithmetic_statement"]
        self.assertEqual(statement["required_ratio"], "n/m=8/5")
        self.assertEqual(statement["diophantine_equation"], "5*n=8*m")
        self.assertEqual(statement["gcd_8_5"], 1)
        self.assertEqual(statement["general_solution"], "m=5*t, n=8*t for positive integer t")
        self.assertEqual(statement["minimal_solution"], {"branch_count_m": 5, "total_binary_exponent_n": 8, "t": 1})
        self.assertFalse(statement["physical_five_branch_derivation_claimed"])

    def test_bounded_scan_matches_multiples(self) -> None:
        scan = self.payload["bounded_exact_match_scan"]
        self.assertEqual(scan["max_branches"], 40)
        self.assertEqual(scan["max_total_binary_exponent"], 80)
        self.assertEqual(scan["exact_match_count"], 8)
        self.assertTrue(scan["all_branch_counts_are_multiples_of_5"])
        self.assertTrue(scan["all_total_binary_exponents_are_multiples_of_8"])
        self.assertEqual(scan["minimal_exact_match_row"]["branch_count_m"], 5)
        self.assertEqual(scan["minimal_exact_match_row"]["total_binary_exponent_n"], 8)
        self.assertLess(abs(scan["minimal_exact_match_row"]["eta_residual_vs_9_5"]), 1e-12)

    def test_minimal_replay_still_needs_selector(self) -> None:
        replay = self.payload["minimal_branch_model_replay"]
        self.assertEqual(replay["minimal_branch_count"], 5)
        self.assertEqual(replay["minimal_total_binary_exponent"], 8)
        self.assertEqual(replay["minimal_labelled_positive_compositions_count"], 35)
        self.assertEqual(replay["minimal_canonical_positive_ledgers"], [[4, 1, 1, 1, 1], [3, 2, 1, 1, 1], [2, 2, 2, 1, 1]])
        self.assertTrue(replay["balanced_ledger_requires_selector_after_minimal_count"])
        self.assertTrue(self.payload["selector_context_replay"]["selector_still_needed_after_minimal_branch_count"])

    def test_candidate_supported_but_model_not_derived(self) -> None:
        candidate = self.payload["candidate_interpretation"]
        self.assertTrue(candidate["supported_by_this_probe"])
        self.assertIn("multiples of five", candidate["content"])
        self.assertIn("Diophantine equation 5n=8m", candidate["why_this_is_more_proof_like"])
        self.assertIn("does not derive", candidate["why_this_is_not_enough"])
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
