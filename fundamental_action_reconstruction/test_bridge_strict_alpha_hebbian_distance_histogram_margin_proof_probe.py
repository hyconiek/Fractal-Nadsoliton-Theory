import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_distance_histogram_margin_proof_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_distance_histogram_margin_proof_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_distance_histogram_margin_proof_report.md"


class StrictAlphaHebbianDistanceHistogramMarginProofProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_HEBBIAN_DISTANCE_HISTOGRAM_MARGIN_PROOF_PROBE__NOT_A_THEOREM",
        )
        self.assertEqual(
            payload["status"],
            "d5-energy-margin-reduced-to-distance-histogram-proof-not-origin-theorem",
        )
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["active_count"], 5)
        self.assertEqual(model["support_count"], 792)
        self.assertEqual(model["histogram_class_count"], 35)
        self.assertEqual(model["teacher_step"], 5)
        self.assertEqual(model["teacher_orbit_size"], 12)
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")

    def test_histogram_proof_core(self):
        core = self.payload["histogram_proof_core"]
        self.assertEqual(core["d5_histogram_d1_to_d6"], [0, 3, 2, 1, 4, 0])
        self.assertEqual(core["d5_histogram_support_count"], 12)
        self.assertTrue(core["d5_histogram_supports_equal_teacher_orbit"])
        self.assertEqual(core["closest_non_d5_histogram_d1_to_d6"], [1, 3, 2, 1, 3, 0])
        self.assertEqual(core["closest_non_d5_histogram_support_count"], 24)
        self.assertEqual(core["closest_non_d5_examples"][0], [0, 1, 3, 5, 10])
        self.assertEqual(core["support_enumeration_reduced_from"], 792)
        self.assertEqual(core["support_enumeration_reduced_to_histogram_classes"], 35)

    def test_histogram_margin_summary(self):
        summary = self.payload["histogram_margin_summary"]
        self.assertEqual(summary["tested_variant_count"], 6)
        self.assertTrue(summary["all_variants_have_positive_histogram_margin"])
        self.assertTrue(summary["all_variants_have_required_stabilizer"])
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
        self.assertEqual(summary["diagonal_shift_by_variant"]["binary_with_diagonal"], "25")
        self.assertEqual(summary["diagonal_shift_by_variant"]["centered_with_diagonal"], "175/12")
        self.assertEqual(summary["diagonal_shift_by_variant"]["bipolar_with_diagonal"], "60")
        self.assertEqual(summary["diagonal_shift_by_variant"]["centered_zero_self"], "0")

    def test_variant_histogram_certificate(self):
        variants = {row["variant_name"]: row for row in self.payload["variant_histogram_certificates"]}
        centered = variants["centered_zero_self"]
        self.assertEqual(centered["unit_stabilizer"], [1, 11])
        self.assertEqual(centered["top_histogram"]["distance_histogram_d1_to_d6"], [0, 3, 2, 1, 4, 0])
        self.assertEqual(centered["top_histogram"]["total_energy"], "55/3")
        self.assertEqual(centered["second_histogram"]["distance_histogram_d1_to_d6"], [1, 3, 2, 1, 3, 0])
        self.assertEqual(centered["second_histogram"]["total_energy"], "31/3")
        self.assertEqual(centered["histogram_gap"], "8")
        self.assertTrue(centered["positive_histogram_margin_certificate"])

        binary_diag = variants["binary_with_diagonal"]
        self.assertEqual(binary_diag["diagonal_shift_for_five_active_supports"], "25")
        self.assertEqual(binary_diag["top_histogram"]["total_energy"], "85")
        self.assertEqual(binary_diag["second_histogram"]["total_energy"], "77")

    def test_interpretation_and_guardrails(self):
        payload = self.payload
        interpretation = payload["candidate_source_interpretation"]
        self.assertIn("compressed from 792 support energies to 35", interpretation["finite_gain"])
        self.assertIn("unique top histogram", interpretation["conditional_positive_result"])
        self.assertIn("does not derive the d5 teacher trace", interpretation["honest_limit"])
        self.assertIn("diagonal-shift replay", interpretation["relation_to_previous_probe"])

        ontology = payload["ontology_guardrail"]
        self.assertIn("nadsoliton itself", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives the d5 teacher/self-record trace", hard_limits)
        self.assertIn("No theorem derives a Hebbian learning law", hard_limits)
        self.assertIn("No theorem derives the pair-distance histogram representation", hard_limits)
        self.assertIn("not an exhaustive theorem", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
