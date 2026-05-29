import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_energy_landscape_d5_maximum_proof_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_energy_landscape_d5_maximum_proof_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_energy_landscape_d5_maximum_proof_report.md"


class StrictAlphaHebbianEnergyLandscapeD5MaximumProofProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_HEBBIAN_ENERGY_LANDSCAPE_D5_MAXIMUM_PROOF_PROBE__NOT_A_THEOREM",
        )
        self.assertEqual(payload["status"], "finite-conditional-d5-global-maximum-proof-not-strict-source")
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["active_count"], 5)
        self.assertEqual(model["support_count"], 792)
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")

    def test_exact_weight_and_energy_landscape(self):
        payload = self.payload
        weight = payload["exact_weight_replay"]
        self.assertEqual(
            weight["first_row_by_cyclic_distance"],
            ["0", "-25/12", "11/12", "-1/12", "-13/12", "23/12", "-25/12", "23/12", "-13/12", "-1/12", "11/12", "-25/12"],
        )
        self.assertTrue(weight["circulant_check"])
        self.assertTrue(weight["diagonal_zero_check"])

        landscape = payload["finite_energy_landscape_certificate"]
        self.assertEqual(landscape["support_count"], 792)
        self.assertEqual(landscape["energy_level_count"], 12)
        self.assertEqual(landscape["maximum_energy"], "55/3")
        self.assertEqual(landscape["maximizer_count"], 12)
        self.assertTrue(landscape["maximizers_equal_d5_teacher_orbit"])
        self.assertTrue(landscape["all_d5_teacher_patterns_are_global_maxima"])
        self.assertEqual(landscape["non_d5_global_maximizer_count"], 0)
        self.assertEqual(landscape["energy_level_histogram"]["55/3"], 12)
        self.assertEqual(landscape["energy_level_histogram"]["-41/3"], 12)

    def test_distance_histogram_proof(self):
        proof = self.payload["distance_histogram_proof_certificate"]
        self.assertEqual(
            proof["distance_weights_d1_to_d6"],
            {"1": "-25/12", "2": "11/12", "3": "-1/12", "4": "-13/12", "5": "23/12", "6": "-25/12"},
        )
        self.assertEqual(proof["histogram_class_count"], 35)
        self.assertTrue(proof["top_histogram_unique"])
        self.assertEqual(proof["top_histogram"]["distance_histogram_d1_to_d6"], [0, 3, 2, 1, 4, 0])
        self.assertEqual(proof["top_histogram"]["support_count"], 12)
        self.assertEqual(proof["top_histogram"]["energy"], "55/3")
        self.assertTrue(proof["direct_energy_equals_histogram_replay_for_all_supports"])
        self.assertEqual(proof["replay_failure_count"], 0)

    def test_energy_refined_ascent_replay(self):
        ascent = self.payload["energy_refined_kwta_ascent_replay"]
        self.assertTrue(ascent["strict_positive_ascent_for_every_non_d5_retained_branch"])
        self.assertEqual(ascent["minimum_positive_delta"], "4")
        self.assertEqual(ascent["delta_histogram"]["0"], 12)
        self.assertEqual(ascent["delta_histogram"]["28"], 132)
        self.assertEqual(ascent["energy_transition_histogram"]["55/3 -> 55/3"], 12)
        self.assertEqual(ascent["violation_count"], 0)

    def test_proof_reading_and_guardrails(self):
        payload = self.payload
        proof = payload["proof_reading"]
        self.assertIn("only global maximizers", proof["finite_theorem"])
        self.assertIn("35 cyclic distance histograms", proof["how_verified"])
        self.assertIn("finite Lyapunov target", proof["relation_to_previous_probe"])
        self.assertIn("still assumed", proof["remaining_gap"])

        ontology = payload["ontology_guardrail"]
        self.assertIn("finite self-recorded resonance patterns", ontology["allowed_reading"])
        self.assertIn("not introduced as a separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives the d5 teacher/self-record trace", hard_limits)
        self.assertIn("finite conditional landscape theorem", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
