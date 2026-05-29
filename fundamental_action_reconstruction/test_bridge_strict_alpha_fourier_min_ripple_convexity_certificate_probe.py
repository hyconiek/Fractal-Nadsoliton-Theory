#!/usr/bin/env python3
"""Tests for the strict-alpha Fourier minimum-ripple convexity certificate probe."""
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_fourier_min_ripple_convexity_certificate_probe.py"
REPORT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_fourier_min_ripple_convexity_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_fourier_min_ripple_convexity_certificate_report.md"


class BridgeStrictAlphaFourierMinRippleConvexityCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(REPORT.read_text(encoding="utf-8"))

    def test_result_kind_and_guardrails(self) -> None:
        self.assertEqual(
            self.payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_FOURIER_MIN_RIPPLE_CONVEXITY_CERTIFICATE_PROBE__NOT_A_THEOREM",
        )
        self.assertIn("No identity K_legacy_ont == K_strict_gate is exported.", self.payload["hard_limits"])
        self.assertIn("No theorem derives the minimum-ripple/coherence-compression principle from strict nadsoliton geometry.", self.payload["hard_limits"])
        self.assertIn("This certificate is conditional on adopting the minimum Fourier ripple selector; it is not a strict-core selector source.", self.payload["hard_limits"])
        self.assertIn("No QW-2191 selector discharge and no ToE closure are claimed.", self.payload["hard_limits"])

    def test_conditional_theorem_statement(self) -> None:
        theorem = self.payload["conditional_theorem_statement"]
        self.assertEqual(theorem["parseval_identity"], "P_N(e)=N*sum_i(e_i^2)-n^2")
        self.assertEqual(theorem["pairwise_smoothing_rule"], "if a>=b+2, replace (a,b) by (a-1,b+1)")
        self.assertEqual(theorem["exact_drop"], "P_N(before)-P_N(after)=2*N*(a-b-1)>0")
        self.assertIn("minimum total non-DC Fourier ripple criterion", theorem["premises"])
        self.assertIn("integer-convexity certificate only", theorem["physical_status"])

    def test_target_certificate_selects_balanced_ledger(self) -> None:
        target = self.payload["target_eta_9_5_certificate"]
        self.assertEqual(target["q_power_product"], "256/243")
        self.assertEqual(target["branch_count"], 5)
        self.assertEqual(target["total_binary_exponent"], 8)
        self.assertLess(abs(target["eta_residual_vs_9_5"]), 1e-12)
        self.assertEqual(target["balanced_ledger"], [2, 2, 2, 1, 1])
        self.assertTrue(target["unique_minimizer"])
        self.assertEqual(target["minimizer"]["ledger"], [2, 2, 2, 1, 1])
        self.assertEqual(target["minimizer"]["sum_squares"], 14)
        self.assertEqual(target["minimizer"]["five_branch_total_non_dc_power"], 6)
        self.assertEqual(target["minimizer"]["z12_total_non_dc_power"], 104)

    def test_smoothing_paths_have_exact_drops(self) -> None:
        paths = self.payload["target_eta_9_5_certificate"]["sample_smoothing_paths"]
        path_41111 = paths["from_4_1_1_1_1"]
        self.assertEqual(len(path_41111), 2)
        self.assertEqual(path_41111[0]["before"], [4, 1, 1, 1, 1])
        self.assertEqual(path_41111[0]["after"], [3, 2, 1, 1, 1])
        self.assertEqual(path_41111[0]["sum_squares_drop"], 4)
        self.assertEqual(path_41111[0]["z12_power_drop"], 48)
        self.assertEqual(path_41111[0]["five_branch_power_drop"], 20)
        self.assertEqual(path_41111[1]["after"], [2, 2, 2, 1, 1])
        self.assertEqual(paths["from_2_2_2_1_1"], [])

    def test_bounded_scan_confirms_balanced_minimizer(self) -> None:
        scan = self.payload["bounded_convexity_scan"]
        self.assertEqual(scan["scan_limits"], {"max_branch_count": 9, "max_total_exponent": 24})
        self.assertEqual(scan["cases_checked"], 180)
        self.assertEqual(scan["canonical_ledgers_checked"], 5583)
        self.assertTrue(scan["all_unique_balanced_minimizers"])
        self.assertEqual(scan["failures"], [])
        target_rows = [row for row in scan["rows"] if row["branch_count"] == 5 and row["total"] == 8]
        self.assertEqual(len(target_rows), 1)
        self.assertEqual(target_rows[0]["balanced_ledger"], [2, 2, 2, 1, 1])
        self.assertEqual(target_rows[0]["min_sum_squares"], 14)
        self.assertTrue(target_rows[0]["unique_minimizer"])

    def test_candidate_supported_but_not_closure(self) -> None:
        candidate = self.payload["candidate_interpretation"]
        self.assertTrue(candidate["supported_by_this_probe"])
        self.assertIn("unique integer convexity minimizer", candidate["content"])
        self.assertIn("pairwise smoothing lemma", candidate["why_this_is_more_proof_like"])
        self.assertIn("does not derive", candidate["why_this_is_not_enough"])
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
