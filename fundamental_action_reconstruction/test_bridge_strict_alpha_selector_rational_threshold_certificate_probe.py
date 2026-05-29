#!/usr/bin/env python3
"""Tests for the strict-alpha rational selector threshold certificate probe."""
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_selector_rational_threshold_certificate_probe.py"
REPORT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_selector_rational_threshold_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_selector_rational_threshold_certificate_report.md"


class BridgeStrictAlphaSelectorRationalThresholdCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(REPORT.read_text(encoding="utf-8"))

    def test_result_kind_and_guardrails(self) -> None:
        self.assertEqual(self.payload["result_kind"], "SCRATCH_STRICT_ALPHA_SELECTOR_RATIONAL_THRESHOLD_CERTIFICATE_PROBE__NOT_A_THEOREM")
        self.assertIn("No identity K_legacy_ont == K_strict_gate is exported.", self.payload["hard_limits"])
        self.assertIn("No legacy physical-role transfer is licensed.", self.payload["hard_limits"])
        self.assertIn("No theorem derives any safe rational orbit-weight exponent p/q with 3^q > 2^(p+q).", self.payload["hard_limits"])
        self.assertIn("No QW-2191 selector discharge and no ToE closure are claimed.", self.payload["hard_limits"])

    def test_exact_rational_certificate_rule(self) -> None:
        rule = self.payload["exact_rational_certificate_rule"]
        self.assertEqual(rule["balanced_wins_iff"], "3^q > 2^(p+q)")
        self.assertEqual(rule["orbit_favoured_wins_iff"], "3^q < 2^(p+q)")
        self.assertTrue(rule["rational_tie_impossible_by_prime_factorization"])
        self.assertAlmostEqual(rule["gamma_c_numeric"], 0.5849625007211562, places=15)

    def test_bounded_scan_best_safe_and_first_fail(self) -> None:
        scan = self.payload["bounded_rational_scan"]
        self.assertEqual(scan["max_denominator"], 64)
        self.assertEqual(scan["tie_count"], 0)
        self.assertEqual(scan["best_safe_below_threshold"]["gamma_label"], "31/53")
        self.assertEqual(scan["first_fail_above_threshold"]["gamma_label"], "24/41")
        self.assertGreater(scan["safe_count"], 0)
        self.assertGreater(scan["fail_count"], 0)

    def test_common_rational_certificates(self) -> None:
        rows = {row["gamma_label"]: row for row in self.payload["common_rational_gamma_certificates"]}
        self.assertEqual(rows["1/2"]["classification"], "SAFE_BALANCED_WINS_BELOW_THRESHOLD")
        self.assertEqual(rows["1/2"]["winner"], [2, 2, 2, 1, 1])
        self.assertEqual(rows["3/5"]["classification"], "FAIL_ORBIT_FAVOURED_WINS_ABOVE_THRESHOLD")
        self.assertEqual(rows["3/5"]["winner"], [3, 2, 1, 1, 1])
        self.assertGreater(rows["31/53"]["integer_certificate_delta_3q_minus_2p_plus_q"], 0)
        self.assertLess(rows["24/41"]["integer_certificate_delta_3q_minus_2p_plus_q"], 0)

    def test_candidate_supported_but_still_conditional(self) -> None:
        candidate = self.payload["candidate_interpretation"]
        self.assertTrue(candidate["supported_by_this_probe"])
        self.assertIn("exact integer certificate", candidate["content"])
        self.assertIn("3^q-2^(p+q)", candidate["why_this_is_more_proof_like"])
        self.assertIn("No strict theorem", candidate["why_this_is_not_enough"])
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
