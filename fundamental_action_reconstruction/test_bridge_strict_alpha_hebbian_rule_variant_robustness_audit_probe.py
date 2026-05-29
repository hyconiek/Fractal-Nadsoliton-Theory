import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_rule_variant_robustness_audit_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_rule_variant_robustness_audit_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_rule_variant_robustness_audit_report.md"


class StrictAlphaHebbianRuleVariantRobustnessAuditProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_HEBBIAN_RULE_VARIANT_ROBUSTNESS_AUDIT_PROBE__NOT_A_THEOREM",
        )
        self.assertEqual(
            payload["status"],
            "d5-hebbian-stabilizer-and-energy-maxima-robust-across-tested-readouts-not-origin-theorem",
        )
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["active_count"], 5)
        self.assertEqual(model["support_count"], 792)
        self.assertEqual(model["teacher_step"], 5)
        self.assertEqual(model["teacher_orbit_size"], 12)
        self.assertEqual(model["minimal_required_subgroup_from_previous_probe"], [1, 11])
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")

    def test_all_variants_pass_robustness_gate(self):
        summary = self.payload["robustness_summary"]
        expected = [
            "binary_with_diagonal",
            "binary_zero_self",
            "centered_with_diagonal",
            "centered_zero_self",
            "bipolar_with_diagonal",
            "bipolar_zero_self",
        ]
        self.assertEqual(summary["tested_variant_count"], 6)
        self.assertEqual(summary["passing_variant_count"], 6)
        self.assertEqual(summary["failing_variant_count"], 0)
        self.assertEqual(summary["passing_variants"], expected)
        self.assertEqual(summary["failing_variants"], [])
        self.assertTrue(summary["all_tested_variants_have_required_stabilizer"])
        self.assertTrue(summary["all_tested_variants_have_unique_d5_global_maxima"])
        self.assertEqual(
            summary["maximum_energy_by_variant"],
            {
                "binary_with_diagonal": "85",
                "binary_zero_self": "60",
                "centered_with_diagonal": "395/12",
                "centered_zero_self": "55/3",
                "bipolar_with_diagonal": "140",
                "bipolar_zero_self": "80",
            },
        )

    def test_variant_certificates(self):
        variants = {row["variant_name"]: row for row in self.payload["variant_certificates"]}
        centered = variants["centered_zero_self"]
        self.assertEqual(centered["unit_stabilizer"], [1, 11])
        self.assertTrue(centered["stabilizer_equals_required_subgroup"])
        self.assertTrue(centered["diagonal_zero_check"])
        self.assertEqual(centered["first_row_by_folded_distance_1_to_6"]["1"], "-25/12")
        self.assertEqual(centered["first_row_by_folded_distance_1_to_6"]["5"], "23/12")
        self.assertEqual(centered["energy_landscape"]["maximum_energy"], "55/3")
        self.assertEqual(centered["energy_landscape"]["maximizer_count"], 12)
        self.assertTrue(centered["energy_landscape"]["maximizers_equal_d5_teacher_orbit"])
        self.assertEqual(centered["energy_landscape"]["non_d5_global_maximizer_count"], 0)

        binary = variants["binary_zero_self"]
        self.assertEqual(binary["unit_stabilizer"], [1, 11])
        self.assertEqual(binary["first_row_by_folded_distance_1_to_6"], {"1": "0", "2": "3", "3": "2", "4": "1", "5": "4", "6": "0"})
        self.assertEqual(binary["energy_landscape"]["maximum_energy"], "60")
        self.assertTrue(binary["passes_d5_robustness_gate"])

    def test_interpretation_and_guardrails(self):
        payload = self.payload
        interpretation = payload["candidate_source_interpretation"]
        self.assertIn("not merely a centered zero-self convention artifact", interpretation["finite_gain"])
        self.assertIn("Hebbian-family associative update", interpretation["conditional_positive_result"])
        self.assertIn("does not derive the teacher trace", interpretation["honest_limit"])
        self.assertIn("nearby Hebbian readout conventions", interpretation["relation_to_previous_probe"])

        ontology = payload["ontology_guardrail"]
        self.assertIn("nadsoliton itself", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives the d5 teacher/self-record trace", hard_limits)
        self.assertIn("No theorem derives a Hebbian learning law", hard_limits)
        self.assertIn("not an exhaustive theorem", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
