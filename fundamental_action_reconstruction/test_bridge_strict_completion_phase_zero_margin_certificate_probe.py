from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_zero_margin_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_zero_margin_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_zero_margin_certificate_report.md"


class StrictCompletionPhaseZeroMarginCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_inputs_and_margin_values(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_PHASE_ZERO_MARGIN_CERTIFICATE__RATIONAL_PARAMETER_ROBUSTNESS",
        )
        self.assertIn("positive-rational-parameter-margin", payload["status"])
        inputs = payload["rational_inputs"]
        self.assertEqual(inputs["pi_lower_bound"]["text"], "333/106")
        self.assertEqual(inputs["pi_upper_bound"]["text"], "355/113")
        self.assertEqual(inputs["strict_omega"]["text"], "743/4000")
        self.assertEqual(inputs["strict_phi"]["text"], "13/80")
        values = payload["margin_values"]
        self.assertEqual(values["left_margin"]["text"], "22897/212000")
        self.assertEqual(values["right_margin"]["text"], "17561/226000")
        self.assertEqual(values["limiting_symmetric_parameter_epsilon"]["text"], "17561/2034000")
        self.assertEqual(values["symmetric_parameter_epsilon"]["text"], "17561/4068000")
        self.assertGreater(values["symmetric_parameter_epsilon"]["decimal"], 0.004)

    def test_worst_case_rows_and_summary(self):
        rows = self.payload["worst_case_inequality_rows"]
        self.assertEqual(len(rows), 4)
        self.assertTrue(all(row["strict_inequality_holds"] for row in rows))
        self.assertTrue(all(row["slack"]["decimal"] >= 0.0 for row in rows))
        by_name = {row["name"]: row for row in rows}
        self.assertGreater(by_name["k0_right_edge_lower_worst_case"]["slack"]["decimal"], 0.0)
        self.assertTrue(by_name["k0_left_edge_upper_worst_case"]["strict_inequality_holds"])
        self.assertTrue(by_name["k1_stays_right_of_d11_worst_case"]["strict_inequality_holds"])

        summary = self.payload["robustness_summary"]
        self.assertTrue(summary["all_margins_positive"])
        self.assertTrue(summary["all_worst_case_inequalities_hold_at_epsilon"])
        self.assertEqual(summary["active_margin_source"], "k0_right_margin/9")
        self.assertTrue(summary["certified_epsilon_is_half_of_limiting_epsilon"])
        self.assertEqual(summary["certified_phase_sign_flip_edges_preserved"], ["1->2", "5->6", "7->8", "9->10"])
        self.assertEqual(summary["certified_phase_sign_pattern_preserved"], [1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1])
        self.assertTrue(summary["matches_phase_zero_report"])

    def test_blocker_context_proof_guardrails_and_markdown(self):
        blockers = self.payload["blocker_context"]
        self.assertIn("rational-pi-intervals", blockers["rational_zero_status"])
        self.assertIn("phase-transport-sign-flips", blockers["phase_zero_status"])
        self.assertIn("strict_phase_frequency_derivation_from_nadsoliton_dynamics", blockers["still_open"])
        self.assertIn("orientation_chi11_source", blockers["still_open"])
        self.assertIn("role_transfer_theorem", blockers["still_open"])

        proof = self.payload["proof_certificate"]
        self.assertIn("7*omega+phi < pi/2", proof["edge_condition_step"])
        self.assertIn("positive rational margins", proof["margin_step"])
        self.assertIn("symmetric epsilon", proof["epsilon_step"])
        self.assertIn("not another zero-location", proof["nonduplication"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("K_strict_gate remains the current live/full", hard_limits)
        self.assertIn("No unqualified identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No proof derives strict omega/phi", hard_limits)
        self.assertIn("No beta_tors -> chi_11 theorem", hard_limits)
        self.assertIn("No QW-2191 selector discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
