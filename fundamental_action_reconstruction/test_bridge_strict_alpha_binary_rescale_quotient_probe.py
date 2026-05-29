#!/usr/bin/env python3
"""Tests for the strict-alpha binary-rescale quotient probe."""
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_binary_rescale_quotient_probe.py"
REPORT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_binary_rescale_quotient_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_binary_rescale_quotient_report.md"


class BridgeStrictAlphaBinaryRescaleQuotientProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(REPORT.read_text(encoding="utf-8"))

    def test_result_kind_and_guardrails(self) -> None:
        self.assertEqual(self.payload["result_kind"], "SCRATCH_STRICT_ALPHA_BINARY_RESCALE_QUOTIENT_PROBE__NOT_A_THEOREM")
        self.assertIn("No identity K_legacy_ont == K_strict_gate is exported.", self.payload["hard_limits"])
        self.assertIn("No legacy physical-role transfer is licensed.", self.payload["hard_limits"])
        self.assertIn("No theorem derives binary-rescale quotienting as a strict nadsoliton gauge principle.", self.payload["hard_limits"])
        self.assertIn("No QW-2191 selector discharge and no ToE closure are claimed.", self.payload["hard_limits"])

    def test_equivalence_rule_is_conditional(self) -> None:
        rule = self.payload["equivalence_rule"]
        self.assertEqual(rule["model_invariant"], "q_model^m = 2^n/D^m")
        self.assertEqual(rule["reduction_step"], "if D even and n>=m, (D,m,n)->(D/2,m,n-m)")
        self.assertEqual(rule["canonical_representative_condition"], "D odd")
        self.assertFalse(rule["physical_gauge_quotient_derived_here"])

    def test_bounded_scan_collapses_to_primitive(self) -> None:
        scan = self.payload["bounded_family_quotient_scan"]
        self.assertEqual(scan["row_count"], 56)
        self.assertEqual(scan["primitive_row_count"], 8)
        self.assertEqual(scan["nonprimitive_row_count"], 48)
        self.assertTrue(scan["all_invariants_preserved"])
        self.assertTrue(scan["all_final_denominators_odd"])
        self.assertTrue(scan["all_reduction_steps_equal_k"])
        self.assertTrue(scan["all_finals_match_expected_primitive"])
        self.assertTrue(scan["all_eta_residuals_exact"])

    def test_sample_k_paths_have_expected_steps(self) -> None:
        samples = {row["binary_rescale_k"]: row for row in self.payload["bounded_family_quotient_scan"]["sample_t1_paths_by_k"]}
        self.assertEqual(samples[0]["canonicalization"]["final"], {"D": 3, "m": 5, "n": 8})
        self.assertEqual(samples[0]["canonicalization"]["step_count"], 0)
        self.assertEqual(samples[1]["canonicalization"]["initial"], {"D": 6, "m": 5, "n": 13})
        self.assertEqual(samples[1]["canonicalization"]["final"], {"D": 3, "m": 5, "n": 8})
        self.assertEqual(samples[1]["canonicalization"]["step_count"], 1)
        self.assertEqual(samples[2]["canonicalization"]["initial"], {"D": 12, "m": 5, "n": 18})
        self.assertEqual(samples[2]["canonicalization"]["final"], {"D": 3, "m": 5, "n": 8})
        self.assertEqual(samples[2]["canonicalization"]["step_count"], 2)

    def test_candidate_supported_but_gauge_not_derived(self) -> None:
        candidate = self.payload["candidate_interpretation"]
        self.assertTrue(candidate["supported_by_this_probe"])
        self.assertIn("binary-rescale equivalent", candidate["content"])
        self.assertIn("invariant-preserving reduction algorithm", candidate["why_this_is_more_proof_like"])
        self.assertIn("does not derive", candidate["why_this_is_not_enough"])
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
