import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_arrow_action_positive_lambda_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_arrow_action_positive_lambda_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_arrow_action_positive_lambda_certificate_report.md"


class StrictAlphaArrowActionPositiveLambdaCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_target_replay(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_ARROW_ACTION_POSITIVE_LAMBDA_CERTIFICATE_PROBE__NOT_A_THEOREM",
        )
        replay = payload["target_identity_replay"]
        self.assertEqual(replay["q_power_product"], "256/243")
        self.assertAlmostEqual(replay["eta_residual_vs_9_5"], 0.0, places=14)
        self.assertEqual(replay["branch_count_m"], 5)
        self.assertEqual(replay["binary_exponent_n"], 8)
        self.assertEqual(replay["balanced_ledger"], [2, 2, 2, 1, 1])

    def test_previous_arrow_replay_and_action_definition(self):
        payload = self.payload
        replay = payload["previous_arrow_action_replay"]
        self.assertIn("self_recorded_arrow_action_lexicographic_selector_report.json", replay["source_report"])
        self.assertEqual(replay["lexicographic_winner"]["winners"], [[2, 2, 2, 1, 1]])
        self.assertTrue(replay["dominance_certificate"]["all_equal_ripple_competitors_have_positive_arrow_penalty"])
        action = payload["scalar_action_definition"]
        self.assertIn("S_lambda(e)=R(e)+lambda*A(e)", action["S_lambda"])
        self.assertIn("lambda>0", action["finite_claim"])
        self.assertIn("max(0, e_{i+1}-e_i)", action["A"])

    def test_zero_and_positive_lambda_scans(self):
        cert = self.payload["finite_positive_lambda_certificate"]
        self.assertEqual(cert["total_positive_ordered_ledgers"], 35)
        self.assertEqual(cert["target_row"]["ordered_ledger"], [2, 2, 2, 1, 1])
        self.assertEqual(cert["target_row"]["ripple_sum_squares_R"], 14)
        self.assertEqual(cert["target_row"]["arrow_increase_penalty_A"], 0)

        zero_scan = cert["lambda_zero_scan"]
        self.assertEqual(zero_scan["lambda"], "0")
        self.assertEqual(zero_scan["winner_count"], 10)
        self.assertFalse(zero_scan["uniquely_selects_balanced"])
        self.assertIn([2, 2, 2, 1, 1], zero_scan["winners"])

        self.assertTrue(cert["all_scanned_positive_lambdas_uniquely_select_balanced"])
        self.assertEqual([scan["lambda"] for scan in cert["positive_lambda_scan"]], ["1/100", "1/10", "1/2", "1", "2", "10", "100"])
        for scan in cert["positive_lambda_scan"]:
            self.assertEqual(scan["winner_count"], 1)
            self.assertEqual(scan["winners"], [[2, 2, 2, 1, 1]])
            self.assertTrue(scan["uniquely_selects_balanced"])

    def test_symbolic_certificate_and_guardrails(self):
        payload = self.payload
        symbolic = payload["finite_positive_lambda_certificate"]["symbolic_certificate"]
        self.assertEqual(symbolic["competitor_count"], 34)
        self.assertEqual(symbolic["equal_ripple_competitor_count"], 9)
        self.assertEqual(symbolic["higher_ripple_competitor_count"], 25)
        self.assertEqual(symbolic["lower_ripple_competitor_count"], 0)
        self.assertEqual(symbolic["equal_ripple_delta_A_values"], [1, 2])
        self.assertEqual(symbolic["higher_ripple_delta_R_values"], [2, 6])
        self.assertEqual(symbolic["higher_ripple_zero_arrow_competitor_count"], 2)
        self.assertEqual(symbolic["min_positive_equal_ripple_delta_A"], 1)
        self.assertEqual(symbolic["min_higher_ripple_delta_R"], 2)
        self.assertTrue(symbolic["all_competitors_satisfy_delta_R_nonnegative"])
        self.assertTrue(symbolic["all_equal_ripple_competitors_have_positive_delta_A"])
        self.assertEqual(symbolic["nonpositive_arrow_equal_ripple_competitors"], [])
        self.assertIn("delta_R + lambda*delta_A", symbolic["positive_lambda_gap_formula"])

        self.assertEqual(
            payload["candidate_interpretation"]["status"],
            "candidate-supported-but-positive-lambda-arrow-action-not-derived",
        )
        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("nadsoliton is treated as information", hard_limits)
        self.assertIn("No theorem derives the arrow-increase penalty", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
