import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_fifth_mode_symbolic_histogram_selector_proof_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_fifth_mode_symbolic_histogram_selector_proof_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_fifth_mode_symbolic_histogram_selector_proof_report.md"


class StrictAlphaHebbianFifthModeSymbolicHistogramSelectorProofProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_HEBBIAN_FIFTH_MODE_SYMBOLIC_HISTOGRAM_SELECTOR_PROOF_PROBE__NOT_A_THEOREM",
        )
        self.assertEqual(payload["status"], "fifth-mode-symbolically-selects-d5-conditionally-not-origin-theorem")
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["active_count"], 5)
        self.assertEqual(model["support_count"], 792)
        self.assertEqual(model["histogram_class_count"], 35)
        self.assertEqual(model["target_mode"], 5)
        self.assertEqual(model["teacher_step"], 5)
        self.assertEqual(model["teacher_orbit_size"], 12)
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")

    def test_symbolic_power_formula(self):
        formula = self.payload["symbolic_power_formula"]
        self.assertEqual(formula["formula"], "P_5(h)=5+2*sum_d h_d*cos(5*pi*d/6)")
        self.assertEqual(formula["coherent_power_field"], "Q(sqrt(3))")
        self.assertEqual(formula["fifth_mode_channel_stabilizer"], [1, 11])
        self.assertIn("no floating-point", formula["comparison_method"])
        self.assertEqual(formula["cos_coefficients_d1_to_d6_as_rational_plus_sqrt3_coeff"][0], {"rational_part": "0", "sqrt3_coefficient": "-1/2"})
        self.assertEqual(formula["cos_coefficients_d1_to_d6_as_rational_plus_sqrt3_coeff"][4], {"rational_part": "0", "sqrt3_coefficient": "1/2"})

    def test_selector_certificate(self):
        certificate = self.payload["selector_certificate"]
        self.assertTrue(certificate["top_histogram_is_d5"])
        self.assertTrue(certificate["second_histogram_is_closest_non_d5"])
        self.assertEqual(certificate["top_histogram"]["distance_histogram_d1_to_d6"], [0, 3, 2, 1, 4, 0])
        self.assertEqual(certificate["top_histogram"]["support_count"], 12)
        self.assertEqual(certificate["top_histogram"]["p5_power_exact"], {"rational_part": "7", "sqrt3_coefficient": "4", "expression": "7 + 4*sqrt3"})
        self.assertEqual(certificate["second_histogram"]["distance_histogram_d1_to_d6"], [1, 3, 2, 1, 3, 0])
        self.assertEqual(certificate["second_histogram"]["p5_power_exact"], {"rational_part": "7", "sqrt3_coefficient": "2", "expression": "7 + 2*sqrt3"})
        self.assertEqual(certificate["d5_histogram_support_count"], 12)
        self.assertTrue(certificate["d5_histogram_supports_equal_teacher_orbit"])
        self.assertEqual(certificate["closest_non_d5_support_count"], 24)
        self.assertEqual(certificate["symbolic_gap_top_minus_second"], {"rational_part": "0", "sqrt3_coefficient": "2", "expression": "2*sqrt3"})
        self.assertTrue(certificate["positive_symbolic_gap"])

    def test_interpretation_and_guardrails(self):
        payload = self.payload
        relation = payload["relation_to_previous_core"]
        self.assertIn("h_5 uniquely selects d5", relation["distance5_pair_count_core"])
        self.assertIn("additional h_1,h_2,h_4,h_6 terms", relation["fifth_mode_readout"])
        self.assertIn("symbolic fifth-resonance", relation["finite_gain"])

        interpretation = payload["candidate_source_interpretation"]
        self.assertIn("Q(sqrt(3))", interpretation["finite_gain"])
        self.assertIn("unique maximizer", interpretation["conditional_positive_result"])
        self.assertIn("does not derive why strict geometry must choose k=5", interpretation["honest_limit"])

        ontology = payload["ontology_guardrail"]
        self.assertIn("nadsoliton itself", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives fifth-mode channel selection", hard_limits)
        self.assertIn("conditional on supplying k=5", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
