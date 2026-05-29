import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_symbolic_mode_competition_no_go_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_symbolic_mode_competition_no_go_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_symbolic_mode_competition_no_go_report.md"


class StrictAlphaHebbianSymbolicModeCompetitionNoGoProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_HEBBIAN_SYMBOLIC_MODE_COMPETITION_NO_GO_PROBE__NOT_A_THEOREM",
        )
        self.assertEqual(payload["status"], "unconstrained-highest-resonance-selects-nyquist-not-d5")
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["active_count"], 5)
        self.assertEqual(model["support_count"], 792)
        self.assertEqual(model["histogram_class_count"], 35)
        self.assertEqual(model["modes_scanned"], [1, 2, 3, 4, 5, 6])
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")

    def test_per_mode_certificates(self):
        modes = self.payload["mode_certificates"]
        self.assertTrue(modes["5"]["top_is_d5"])
        self.assertEqual(modes["5"]["top_histogram"]["power_exact"], {"rational_part": "7", "sqrt3_coefficient": "4", "expression": "7 + 4*sqrt3"})
        self.assertEqual(modes["5"]["top_minus_second_gap_exact"], {"rational_part": "0", "sqrt3_coefficient": "2", "expression": "2*sqrt3"})
        self.assertTrue(modes["1"]["top_is_contiguous"])
        self.assertEqual(modes["1"]["top_histogram"]["distance_histogram_d1_to_d6"], [4, 3, 2, 1, 0, 0])
        self.assertTrue(modes["6"]["top_is_nyquist"])
        self.assertEqual(modes["6"]["top_histogram"]["power_exact"], {"rational_part": "25", "sqrt3_coefficient": "0", "expression": "25"})
        self.assertEqual(modes["6"]["channel_stabilizer"], [1, 5, 7, 11])

    def test_mode_competition_no_go(self):
        certificate = self.payload["mode_competition_certificate"]
        self.assertEqual(certificate["global_winner_mode"], 6)
        self.assertEqual(certificate["global_winner_power_exact"], {"rational_part": "25", "sqrt3_coefficient": "0", "expression": "25"})
        self.assertTrue(certificate["global_winner_is_nyquist"])
        self.assertEqual(certificate["fifth_mode_d5_power_exact"], {"rational_part": "7", "sqrt3_coefficient": "4", "expression": "7 + 4*sqrt3"})
        self.assertEqual(certificate["global_winner_minus_fifth_mode_d5_exact"], {"rational_part": "18", "sqrt3_coefficient": "-4", "expression": "18 - 4*sqrt3"})
        self.assertTrue(certificate["global_winner_strictly_exceeds_fifth_mode_d5"])
        self.assertTrue(certificate["fifth_mode_selects_d5_conditionally"])
        self.assertFalse(certificate["unconstrained_highest_resonance_selects_d5"])

    def test_interpretation_and_guardrails(self):
        payload = self.payload
        interpretation = payload["candidate_source_interpretation"]
        self.assertIn("conditional fifth-mode d5 selector", interpretation["finite_gain"])
        self.assertIn("k=5 is supplied", interpretation["conditional_positive_result"])
        self.assertIn("k=6/Nyquist wins", interpretation["negative_result"])
        self.assertIn("does not derive an internal source", interpretation["honest_limit"])

        ontology = payload["ontology_guardrail"]
        self.assertIn("nadsoliton itself", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives fifth-mode channel selection", hard_limits)
        self.assertIn("k=6 Nyquist winner", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
