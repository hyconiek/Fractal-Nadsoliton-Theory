import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_energy_margin_robustness_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_energy_margin_robustness_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_energy_margin_robustness_certificate_report.md"


class StrictAlphaHebbianEnergyMarginRobustnessCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_HEBBIAN_ENERGY_MARGIN_ROBUSTNESS_CERTIFICATE_PROBE__NOT_A_THEOREM",
        )
        self.assertEqual(
            payload["status"],
            "positive-d5-energy-margin-across-tested-hebbian-readouts-not-origin-theorem",
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

    def test_margin_summary(self):
        summary = self.payload["margin_summary"]
        self.assertEqual(summary["tested_variant_count"], 6)
        self.assertTrue(summary["all_tested_variants_have_positive_margin"])
        self.assertTrue(summary["all_tested_variants_have_required_stabilizer"])
        self.assertEqual(summary["minimum_gap"], "8")
        self.assertEqual(
            summary["gap_by_variant"],
            {
                "binary_with_diagonal": "8",
                "binary_zero_self": "8",
                "centered_with_diagonal": "8",
                "centered_zero_self": "8",
                "bipolar_with_diagonal": "32",
                "bipolar_zero_self": "32",
            },
        )
        self.assertTrue(all(value == 24 for value in summary["best_non_d5_count_by_variant"].values()))
        self.assertEqual(summary["support_energy_perturbation_safe_radius_by_variant_strict"]["centered_zero_self"], "4")
        self.assertEqual(summary["entrywise_weight_perturbation_sufficient_bound_by_variant_strict"]["centered_zero_self"], "4/25")
        self.assertEqual(summary["entrywise_weight_perturbation_sufficient_bound_by_variant_strict"]["bipolar_zero_self"], "16/25")

    def test_variant_margin_certificate(self):
        variants = {row["variant_name"]: row for row in self.payload["margin_certificates"]}
        centered = variants["centered_zero_self"]
        self.assertEqual(centered["unit_stabilizer"], [1, 11])
        self.assertEqual(centered["maximum_energy"], "55/3")
        self.assertEqual(centered["best_non_d5_energy"], "31/3")
        self.assertEqual(centered["energy_gap_to_best_non_d5"], "8")
        self.assertEqual(centered["d5_maximizer_count"], 12)
        self.assertTrue(centered["d5_maximizers_equal_teacher_orbit"])
        self.assertEqual(centered["best_non_d5_competitor_count"], 24)
        self.assertEqual(centered["best_non_d5_competitor_examples"][0], [0, 1, 3, 5, 10])
        self.assertEqual(centered["non_d5_gap_histogram"]["8"], 24)
        self.assertTrue(centered["positive_margin_certificate"])

        bipolar = variants["bipolar_zero_self"]
        self.assertEqual(bipolar["maximum_energy"], "80")
        self.assertEqual(bipolar["best_non_d5_energy"], "48")
        self.assertEqual(bipolar["energy_gap_to_best_non_d5"], "32")
        self.assertEqual(bipolar["support_energy_perturbation_safe_radius_strict"], "16")

    def test_interpretation_and_guardrails(self):
        payload = self.payload
        interpretation = payload["candidate_source_interpretation"]
        self.assertIn("exact positive gap", interpretation["finite_gain"])
        self.assertIn("cannot change the d5 winner", interpretation["conditional_positive_result"])
        self.assertIn("does not derive the teacher trace", interpretation["honest_limit"])
        self.assertIn("sufficient perturbation bounds", interpretation["relation_to_previous_probe"])

        ontology = payload["ontology_guardrail"]
        self.assertIn("nadsoliton itself", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives the d5 teacher/self-record trace", hard_limits)
        self.assertIn("No theorem derives a Hebbian learning law", hard_limits)
        self.assertIn("No theorem derives the perturbation model", hard_limits)
        self.assertIn("not an exhaustive theorem", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
