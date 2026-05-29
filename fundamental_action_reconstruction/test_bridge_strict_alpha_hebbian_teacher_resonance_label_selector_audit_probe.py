import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_teacher_resonance_label_selector_audit_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_teacher_resonance_label_selector_audit_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_teacher_resonance_label_selector_audit_report.md"


class StrictAlphaHebbianTeacherResonanceLabelSelectorAuditProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_HEBBIAN_TEACHER_RESONANCE_LABEL_SELECTOR_AUDIT_PROBE__NOT_A_THEOREM",
        )
        self.assertEqual(
            payload["status"],
            "conditional-nonnyquist-unit-resonance-label-selects-d5-not-strict-source",
        )
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["active_count"], 5)
        self.assertEqual(model["support_count"], 792)
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")

    def test_orbit_resonance_labels(self):
        certs = self.payload["orbit_class_certificates"]
        self.assertEqual(certs["contiguous_step_1_11"]["learned_leading_modes"], [1])
        self.assertEqual(certs["parity_minus_one_step_2_10"]["learned_leading_modes"], [6])
        self.assertEqual(certs["fifth_step_d5_step_5_7"]["learned_leading_modes"], [5])
        self.assertAlmostEqual(certs["contiguous_step_1_11"]["learned_leading_score"], 13.928203230276, places=12)
        self.assertAlmostEqual(certs["fifth_step_d5_step_5_7"]["learned_leading_score"], 13.928203230276, places=12)
        self.assertEqual(certs["parity_minus_one_step_2_10"]["learned_leading_score"], 25.0)
        self.assertTrue(certs["contiguous_step_1_11"]["is_unit_step_class"])
        self.assertFalse(certs["parity_minus_one_step_2_10"]["is_unit_step_class"])
        self.assertTrue(certs["fifth_step_d5_step_5_7"]["is_unit_step_class"])

    def test_selector_stack_audit(self):
        stack = self.payload["selector_stack_audit"]
        self.assertEqual(
            stack["hebbian_self_max_only"]["survivors"],
            ["contiguous_step_1_11", "parity_minus_one_step_2_10", "fifth_step_d5_step_5_7"],
        )
        self.assertFalse(stack["hebbian_self_max_only"]["selects_d5_uniquely"])
        self.assertEqual(stack["unit_step_filter"]["survivors"], ["contiguous_step_1_11", "fifth_step_d5_step_5_7"])
        self.assertFalse(stack["unit_step_filter"]["selects_d5_uniquely"])
        self.assertEqual(stack["unqualified_highest_learned_label"]["winners"], ["parity_minus_one_step_2_10"])
        self.assertFalse(stack["unqualified_highest_learned_label"]["selects_d5_uniquely"])
        self.assertIn("k=6 parity/Nyquist", stack["unqualified_highest_learned_label"]["verdict"])
        self.assertEqual(stack["non_nyquist_unit_highest_label"]["winners"], ["fifth_step_d5_step_5_7"])
        self.assertTrue(stack["non_nyquist_unit_highest_label"]["selects_d5_uniquely"])
        self.assertEqual(stack["explicit_fifth_mode_lock"]["winners"], ["fifth_step_d5_step_5_7"])

    def test_proof_reading_and_guardrails(self):
        payload = self.payload
        proof = payload["proof_reading"]
        self.assertIn("k=5 rather than k=1", proof["finite_gain"])
        self.assertIn("k=6 parity/Nyquist", proof["negative_result"])
        self.assertIn("selects d5 uniquely", proof["conditional_positive_result"])
        self.assertIn("not derived", proof["remaining_gap"])

        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives the d5 teacher/self-record trace", hard_limits)
        self.assertIn("No theorem derives the primitive non-Nyquist/unit-channel premise", hard_limits)
        self.assertIn("Unqualified highest resonance remains insufficient", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
