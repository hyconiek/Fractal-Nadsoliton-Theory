import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_highest_resonance_mode_lock_audit_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_highest_resonance_mode_lock_audit_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_highest_resonance_mode_lock_audit_report.md"


class StrictAlphaHighestResonanceModeLockAuditProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_finite_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_HIGHEST_RESONANCE_MODE_LOCK_AUDIT_PROBE__NOT_A_THEOREM",
        )
        self.assertEqual(payload["status"], "candidate-supported-only-with-fifth-mode-lock")
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["support_size"], 5)
        self.assertEqual(model["support_count"], 792)
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")

    def test_fifth_mode_lock_selects_d5_supports(self):
        cert = self.payload["mode_locked_certificate"]
        self.assertEqual(cert["target_mode"], 5)
        self.assertEqual(cert["conjugate_target_mode"], 7)
        self.assertEqual(cert["target_mode_scan"]["winner_count"], 12)
        self.assertEqual(cert["conjugate_target_mode_scan"]["winner_count"], 12)
        self.assertAlmostEqual(cert["target_mode_scan"]["maximum_power"], 13.928203230276, places=12)
        self.assertTrue(cert["target_mode_winners_equal_fifth_step_supports"])
        self.assertTrue(cert["conjugate_mode_winners_equal_fifth_step_supports"])
        self.assertEqual(cert["target_mode_scan"]["winners"], cert["fifth_step_supports"])

    def test_unqualified_highest_resonance_obstruction(self):
        audit = self.payload["unqualified_highest_resonance_audit"]
        global_scan = audit["global_single_mode_scan"]
        self.assertEqual(global_scan["modes_checked"], [1, 2, 3, 4, 5, 6])
        self.assertEqual(global_scan["winning_modes"], [6])
        self.assertEqual(global_scan["global_maximum_power"], 25.0)
        self.assertEqual(global_scan["winner_count"], 12)
        self.assertFalse(global_scan["is_d5_selector"])
        self.assertFalse(audit["unqualified_highest_resonance_selects_d5"])
        self.assertTrue(audit["mode_1_winners_are_contiguous_supports"])
        self.assertTrue(audit["mode_5_winners_are_d5_supports"])
        self.assertIn("k=6 parity/Nyquist", audit["obstruction"])

    def test_coprime_packet_and_guardrails(self):
        payload = self.payload
        packet = payload["coprime_packet_audit"]
        self.assertEqual(packet["packet_scan"]["modes"], [1, 5, 7, 11])
        self.assertEqual(packet["winner_count"], 72)
        self.assertTrue(packet["contains_all_d5_supports"])
        self.assertFalse(packet["is_unique_d5_selector"])

        interpretation = payload["candidate_interpretation"]
        self.assertIn("highest resonance", interpretation["honest_gain"])
        self.assertIn("unqualified highest resonance", interpretation["honest_limit"])
        self.assertIn("primordial information", interpretation["ontology_guardrail"])
        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives the fifth-mode resonance lock", hard_limits)
        self.assertIn("Unqualified highest resonance is explicitly not enough", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
