import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_unit_label_aut_obstruction_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_unit_label_aut_obstruction_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_unit_label_aut_obstruction_report.md"


class StrictAlphaHebbianUnitLabelAutObstructionProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_HEBBIAN_UNIT_LABEL_AUT_OBSTRUCTION_PROBE__NOT_A_THEOREM",
        )
        self.assertEqual(payload["status"], "unit-highest-label-d5-selector-is-aut-breaking-not-strict")
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["active_count"], 5)
        self.assertEqual(model["support_count"], 792)
        self.assertEqual(model["automorphism_units"], [1, 5, 7, 11])
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")

    def test_class_label_replay(self):
        labels = self.payload["class_label_replay"]
        self.assertEqual(labels["contiguous_step_1_11"]["learned_leading_labels"], [1])
        self.assertEqual(labels["parity_minus_one_step_2_10"]["learned_leading_labels"], [6])
        self.assertEqual(labels["fifth_step_d5_step_5_7"]["learned_leading_labels"], [5])
        self.assertEqual(labels["contiguous_step_1_11"]["teacher_support_count"], 12)
        self.assertEqual(labels["fifth_step_d5_step_5_7"]["teacher_support_count"], 12)

    def test_aut_action_and_label_equivariance(self):
        audit = self.payload["aut_z12_action_audit"]
        self.assertTrue(audit["all_class_images_identified"])
        self.assertTrue(audit["all_label_equivariance_checks_pass"])
        self.assertTrue(audit["unit_classes_in_same_aut_orbit"])
        self.assertTrue(audit["parity_class_separate_under_units"])
        rows = audit["rows"]
        map_row = next(row for row in rows if row["source_class"] == "contiguous_step_1_11" and row["unit"] == 5)
        self.assertEqual(map_row["image_class"], "fifth_step_d5_step_5_7")
        self.assertEqual(map_row["source_leading_labels"], [1])
        self.assertEqual(map_row["transformed_leading_labels"], [5])
        self.assertEqual(map_row["image_class_leading_labels"], [5])

    def test_numeric_label_selector_obstruction(self):
        obstruction = self.payload["numeric_label_selector_obstruction"]
        self.assertEqual(obstruction["unit_pair_orbit"], ["contiguous_step_1_11", "fifth_step_d5_step_5_7"])
        self.assertEqual(obstruction["aut_connecting_unit"], 5)
        self.assertEqual(obstruction["contiguous_maps_to"], "fifth_step_d5_step_5_7")
        self.assertEqual(obstruction["label_1_maps_to"], 5)
        self.assertEqual(obstruction["numeric_highest_label_winner"], "fifth_step_d5_step_5_7")
        self.assertFalse(obstruction["is_numeric_highest_label_aut_invariant"])
        self.assertIn("non-invariant numeric label order", obstruction["reason"])

    def test_proof_reading_and_guardrails(self):
        payload = self.payload
        proof = payload["proof_reading"]
        self.assertIn("label-equivariantly", proof["finite_gain"])
        self.assertIn("not Aut-invariant", proof["negative_result"])
        self.assertIn("symmetry-breaking premise", proof["relation_to_previous_probe"])
        self.assertIn("orientation", proof["remaining_gap"])

        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives the d5 teacher/self-record trace", hard_limits)
        self.assertIn("No theorem derives a strict orientation/generator/label-order source", hard_limits)
        self.assertIn("Aut(Z_12)-breaking", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
