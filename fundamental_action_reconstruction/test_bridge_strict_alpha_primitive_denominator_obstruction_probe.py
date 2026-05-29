#!/usr/bin/env python3
"""Tests for the strict-alpha primitive denominator obstruction probe."""
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_primitive_denominator_obstruction_probe.py"
REPORT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_primitive_denominator_obstruction_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_primitive_denominator_obstruction_report.md"


class BridgeStrictAlphaPrimitiveDenominatorObstructionProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(REPORT.read_text(encoding="utf-8"))

    def test_result_kind_and_guardrails(self) -> None:
        self.assertEqual(self.payload["result_kind"], "SCRATCH_STRICT_ALPHA_PRIMITIVE_DENOMINATOR_OBSTRUCTION_PROBE__NOT_A_THEOREM")
        self.assertIn("No identity K_legacy_ont == K_strict_gate is exported.", self.payload["hard_limits"])
        self.assertIn("No legacy physical-role transfer is licensed.", self.payload["hard_limits"])
        self.assertIn("No theorem derives the primitive denominator D=3 from strict nadsoliton geometry.", self.payload["hard_limits"])
        self.assertIn("No QW-2191 selector discharge and no ToE closure are claimed.", self.payload["hard_limits"])

    def test_exact_arithmetic_statement(self) -> None:
        statement = self.payload["exact_arithmetic_statement"]
        self.assertEqual(statement["required_denominator_form"], "D=3*2^(n/m-8/5)")
        self.assertEqual(statement["integer_D_condition"], "n/m-8/5 must be a nonnegative integer k")
        self.assertEqual(statement["general_solution"], "D=3*2^k, m=5*t, n=(8+5*k)*t for integers k>=0,t>=1")
        self.assertEqual(statement["primitive_solution"], "k=0 gives D=3, m=5*t, n=8*t")
        self.assertFalse(statement["denominator_3_derived_from_target_alone"])

    def test_bounded_family_scan_counts_and_exactness(self) -> None:
        scan = self.payload["bounded_family_scan"]
        self.assertEqual(scan["max_binary_rescale_k"], 6)
        self.assertEqual(scan["max_refinement_t"], 8)
        self.assertEqual(scan["exact_family_row_count"], 56)
        self.assertEqual(scan["primitive_row_count"], 8)
        self.assertEqual(scan["nonprimitive_binary_rescale_row_count"], 48)
        self.assertTrue(scan["all_rows_eta_exact"])
        self.assertEqual(scan["non_family_exact_rows_found"], [])

    def test_denominator_summary_contains_rescaled_families(self) -> None:
        summary = {row["binary_rescale_k"]: row for row in self.payload["bounded_family_scan"]["denominator_summary"]}
        self.assertEqual(summary[0]["denominator_D"], 3)
        self.assertEqual(summary[0]["n_over_m"], "8/5")
        self.assertEqual(summary[1]["denominator_D"], 6)
        self.assertEqual(summary[1]["n_over_m"], "13/5")
        self.assertEqual(summary[2]["denominator_D"], 12)
        self.assertEqual(summary[2]["n_over_m"], "18/5")
        self.assertTrue(summary[0]["minimal_t_row"]["primitive_no_binary_rescale"])
        self.assertFalse(summary[1]["minimal_t_row"]["primitive_no_binary_rescale"])

    def test_candidate_supported_but_primitive_not_derived(self) -> None:
        candidate = self.payload["candidate_interpretation"]
        self.assertTrue(candidate["supported_by_this_probe"])
        self.assertIn("D=3*2^k", candidate["content"])
        self.assertIn("primitive denominator selection", candidate["why_this_is_more_proof_like"])
        self.assertIn("does not derive", candidate["why_this_is_not_enough"])
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
