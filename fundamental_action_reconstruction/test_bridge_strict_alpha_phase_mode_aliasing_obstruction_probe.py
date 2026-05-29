import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_phase_mode_aliasing_obstruction_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_phase_mode_aliasing_obstruction_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_phase_mode_aliasing_obstruction_report.md"


class StrictAlphaPhaseModeAliasingObstructionProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_target_replay(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_PHASE_MODE_ALIASING_OBSTRUCTION_PROBE__NOT_A_THEOREM",
        )
        self.assertEqual(payload["target_identity_replay"]["q_power_product"], "256/243")
        self.assertAlmostEqual(payload["target_identity_replay"]["eta_residual_vs_9_5"], 0.0, places=14)
        self.assertEqual(payload["target_identity_replay"]["forward_assignment"], [2, 2, 2, 1, 1])

    def test_previous_phase_certificate_replay_and_theorem_statement(self):
        payload = self.payload
        replay = payload["previous_phase_source_certificate_replay"]
        self.assertTrue(replay["all_rows_recovered_by_all_coprime_modes"])
        self.assertEqual(replay["unique_predicted_orientation_source_pairs_using_mode_1"], 24)
        statement = payload["aliasing_theorem_statement"]
        self.assertIn("k*s modulo 12", statement["phase_observable"])
        self.assertIn("gcd(k,12)=1", statement["unique_source_condition"])
        self.assertIn("gcd(k,12)=g>1", statement["aliasing_condition"])

    def test_gcd_phase_class_counts(self):
        scan = self.payload["phase_mode_aliasing_scan"]
        self.assertEqual(
            scan["expected_phase_class_count_by_mode"],
            {"0": 1, "1": 12, "2": 6, "3": 4, "4": 3, "5": 12, "6": 2, "7": 12, "8": 3, "9": 4, "10": 6, "11": 12},
        )
        self.assertEqual(scan["source_complete_modes"], [1, 5, 7, 11])
        self.assertEqual(scan["coprime_modes_mod_12"], [1, 5, 7, 11])
        self.assertTrue(scan["all_observed_counts_match_gcd_formula"])
        self.assertTrue(scan["all_alias_sizes_match_gcd_formula"])
        summaries = scan["summaries_by_orientation_and_mode"]
        self.assertTrue(summaries["-1"]["1"]["all_sources_uniquely_recovered"])
        self.assertFalse(summaries["-1"]["2"]["all_sources_uniquely_recovered"])
        self.assertEqual(summaries["-1"]["2"]["candidate_set_sizes"], [2])
        self.assertEqual(summaries["1"]["6"]["candidate_set_sizes"], [6])
        self.assertEqual(summaries["1"]["0"]["candidate_set_sizes"], [12])

    def test_guardrails(self):
        payload = self.payload
        self.assertEqual(
            payload["candidate_interpretation"]["status"],
            "candidate-supported-but-phase-origin-still-not-derived",
        )
        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("phase-mode aliasing classification", hard_limits)
        self.assertIn("Non-coprime modes", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
