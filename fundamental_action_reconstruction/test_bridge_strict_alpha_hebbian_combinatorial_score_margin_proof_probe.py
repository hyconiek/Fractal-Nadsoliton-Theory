import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_combinatorial_score_margin_proof_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_combinatorial_score_margin_proof_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_combinatorial_score_margin_proof_report.md"


class StrictAlphaHebbianCombinatorialScoreMarginProofProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_HEBBIAN_COMBINATORIAL_SCORE_MARGIN_PROOF_PROBE__NOT_A_THEOREM",
        )
        self.assertEqual(
            payload["status"],
            "d5-margin-reduced-to-integer-histogram-score-not-origin-theorem",
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

    def test_score_definition(self):
        definition = self.payload["score_definition"]
        self.assertEqual(definition["coefficients_d1_to_d6"], [0, 3, 2, 1, 4, 0])
        self.assertEqual(definition["formula"], "C(h)=3*h_2 + 2*h_3 + h_4 + 4*h_5")
        self.assertEqual(definition["unit_stabilizer"], [1, 11])
        self.assertIn("positive scale", definition["replay_note"])

    def test_score_margin_certificate(self):
        certificate = self.payload["score_margin_certificate"]
        self.assertEqual(certificate["top_histogram"]["distance_histogram_d1_to_d6"], [0, 3, 2, 1, 4, 0])
        self.assertEqual(certificate["top_histogram"]["support_count"], 12)
        self.assertEqual(certificate["top_histogram"]["combinatorial_score"], 30)
        self.assertEqual(certificate["second_histogram"]["distance_histogram_d1_to_d6"], [1, 3, 2, 1, 3, 0])
        self.assertEqual(certificate["second_histogram"]["support_count"], 24)
        self.assertEqual(certificate["second_histogram"]["combinatorial_score"], 26)
        self.assertEqual(certificate["score_gap"], 4)
        self.assertEqual(certificate["d5_histogram_support_count"], 12)
        self.assertTrue(certificate["d5_histogram_supports_equal_teacher_orbit"])
        self.assertEqual(certificate["closest_non_d5_support_count"], 24)
        self.assertEqual(certificate["closest_non_d5_examples"][0], [0, 1, 3, 5, 10])
        self.assertTrue(certificate["top_histogram_unique_by_score"])
        self.assertEqual(certificate["score_level_histogram"]["30"], 1)

    def test_energy_gap_replay(self):
        replay = self.payload["energy_gap_replay"]
        self.assertEqual(replay["energy_gap_scale_by_readout_family"], {"binary": 2, "centered": 2, "bipolar": 8})
        self.assertEqual(
            replay["variant_gap_replay"],
            {
                "binary_with_diagonal": 8,
                "binary_zero_self": 8,
                "centered_with_diagonal": 8,
                "centered_zero_self": 8,
                "bipolar_with_diagonal": 32,
                "bipolar_zero_self": 32,
            },
        )
        self.assertTrue(replay["matches_previous_histogram_margin_probe"])

    def test_interpretation_and_guardrails(self):
        payload = self.payload
        interpretation = payload["candidate_source_interpretation"]
        self.assertIn("integer inequality", interpretation["finite_gain"])
        self.assertIn("positive scale factors", interpretation["conditional_positive_result"])
        self.assertIn("does not derive that setup", interpretation["honest_limit"])
        self.assertIn("underlying integer score certificate", interpretation["relation_to_previous_probe"])

        ontology = payload["ontology_guardrail"]
        self.assertIn("nadsoliton itself", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives the d5 teacher/self-record trace", hard_limits)
        self.assertIn("No theorem derives a Hebbian learning law", hard_limits)
        self.assertIn("No theorem derives the integer histogram score", hard_limits)
        self.assertIn("not exhaustive over all learning laws", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
