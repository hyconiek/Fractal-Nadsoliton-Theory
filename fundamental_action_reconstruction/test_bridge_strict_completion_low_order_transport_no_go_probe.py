from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_low_order_transport_no_go_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_low_order_transport_no_go_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_low_order_transport_no_go_report.md"


class StrictCompletionLowOrderTransportNoGoProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_domain_and_input_summary(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_LOW_ORDER_TRANSPORT_NO_GO__FINITE_Z12_NOT_FULL_DYNAMICAL_NO_GO",
        )
        self.assertIn("low-order-autonomous", payload["status"])
        self.assertEqual(payload["domain"]["node_count"], 12)
        self.assertEqual(payload["domain"]["edge_count"], 11)
        input_summary = payload["transport_input_summary"]
        self.assertEqual(input_summary["negative_edge_count"], 4)
        self.assertEqual(input_summary["negative_edges"], ["1->2", "5->6", "7->8", "9->10"])
        self.assertEqual(input_summary["positive_edge_count"], 7)
        self.assertGreater(input_summary["edge_ratio_spread"], 9.0)

    def test_no_go_checks_fail_as_expected(self):
        checks = self.payload["no_go_checks"]
        self.assertFalse(checks["positive_damping_only_transport"]["exact_pass"])
        self.assertIn("four adjacent ratios", checks["positive_damping_only_transport"]["failure_reason"])

        first_order = checks["constant_first_order_multiplier"]
        self.assertFalse(first_order["exact_pass"])
        self.assertAlmostEqual(first_order["best_r_least_squares"], 0.29307548057349686)
        self.assertGreater(first_order["l2_residual"], 0.6)

        second_order = checks["constant_second_order_linear_recurrence"]
        self.assertFalse(second_order["exact_pass"])
        self.assertAlmostEqual(second_order["best_a_least_squares"], -0.15906840053672286)
        self.assertAlmostEqual(second_order["best_b_least_squares"], -0.07646138093067381)
        self.assertGreater(second_order["l2_residual"], 0.05)

        affine_log = checks["affine_log_envelope"]
        self.assertFalse(affine_log["exact_pass"])
        self.assertGreater(affine_log["l2_residual"], 3.0)

        sign_law = checks["short_period_edge_sign_law"]
        self.assertFalse(sign_law["exact_pass"])
        self.assertEqual([row["period"] for row in sign_law["period_rows"]], [1, 2, 3, 4, 5])
        self.assertTrue(all(not row["is_periodic_on_audited_edges"] for row in sign_law["period_rows"]))

    def test_no_go_summary_blockers_proof_and_guardrails(self):
        summary = self.payload["no_go_summary"]
        self.assertTrue(summary["positive_damping_only_fails"])
        self.assertTrue(summary["constant_first_order_fails"])
        self.assertTrue(summary["constant_second_order_fails"])
        self.assertTrue(summary["affine_log_envelope_fails"])
        self.assertTrue(summary["short_period_edge_sign_law_fails"])
        self.assertIn("nontrivial enough", summary["sharpener"])

        blockers = self.payload["blocker_context"]
        self.assertIn("unique-finite-node-transport", blockers["transport_cocycle_status"])
        self.assertIn("all-three-explicit-factors", blockers["necessity_status"])
        self.assertIn("strict_transport_derivation remains open", blockers["refined_blocker"])
        self.assertIn("orientation_chi11_source", blockers["still_open"])
        self.assertIn("role_transfer_theorem", blockers["still_open"])

        proof = self.payload["proof_certificate"]
        self.assertIn("positive damping-only", proof["positive_damping_no_go"])
        self.assertIn("all edge ratios are equal", proof["first_order_no_go"])
        self.assertIn("nonzero residual", proof["second_order_no_go"])
        self.assertIn("not a theorem against all possible strict dynamics", proof["theoretical_limit"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("No full no-go theorem", hard_limits)
        self.assertIn("No proof derives the completion transport", hard_limits)
        self.assertIn("No beta_tors -> chi_11 theorem", hard_limits)
        self.assertIn("No QW-2191 selector discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
