import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_resonance_learning_audit_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_resonance_learning_audit_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_resonance_learning_audit_report.md"


class StrictAlphaHebbianResonanceLearningAuditProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_finite_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_HEBBIAN_RESONANCE_LEARNING_AUDIT_PROBE__NOT_A_THEOREM",
        )
        self.assertEqual(
            payload["status"],
            "hebbian-amplification-is-conditional-not-a-fifth-mode-source-theorem",
        )
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["support_size"], 5)
        self.assertEqual(model["support_count"], 792)
        self.assertEqual(model["independent_non_dc_modes"], [1, 2, 3, 4, 5, 6])
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")

    def test_dataset_hebbian_certificates(self):
        certs = self.payload["dataset_certificates"]
        all_supports = certs["unbiased_all_supports"]
        self.assertEqual(all_supports["pattern_count"], 792)
        self.assertEqual(all_supports["leading_non_dc_modes"], [1, 2, 3, 4, 5, 6])
        for score in all_supports["average_non_dc_mode_scores"].values():
            self.assertAlmostEqual(score, 3.181818181818, places=12)

        d5 = certs["fifth_step_d5_orbit"]
        self.assertEqual(d5["pattern_count"], 12)
        self.assertEqual(d5["leading_non_dc_modes"], [5])
        self.assertEqual(d5["unique_leading_mode"], 5)
        self.assertAlmostEqual(d5["average_non_dc_mode_scores"]["5"], 13.928203230276, places=12)

        contiguous = certs["contiguous_orbit"]
        self.assertEqual(contiguous["leading_non_dc_modes"], [1])
        self.assertEqual(contiguous["unique_leading_mode"], 1)

        nyquist = certs["parity_nyquist_orbit"]
        self.assertEqual(nyquist["leading_non_dc_modes"], [6])
        self.assertEqual(nyquist["unique_leading_mode"], 6)
        self.assertAlmostEqual(nyquist["average_non_dc_mode_scores"]["6"], 25.0, places=12)

    def test_individual_support_basin_audit(self):
        audit = self.payload["individual_support_basin_audit"]
        self.assertEqual(audit["support_count"], 792)
        self.assertEqual(audit["winner_class_count"], 10)
        self.assertEqual(audit["unique_winner_count"], 696)
        self.assertEqual(audit["tied_winner_count"], 96)
        self.assertEqual(
            audit["unique_winner_counts_by_mode"],
            {"1": 84, "2": 120, "3": 108, "4": 120, "5": 84, "6": 180},
        )
        self.assertEqual(audit["target_mode_unique_count"], 84)
        self.assertEqual(audit["target_mode_including_ties_count"], 144)
        self.assertEqual(audit["target_mode_unique_fraction"], "84/792")
        self.assertEqual(audit["target_mode_including_ties_fraction"], "144/792")
        self.assertEqual(audit["tied_winner_classes"], {"1/2/4/5": 36, "1/3/5": 24, "2/4": 24, "2/6": 12})

    def test_interpretation_and_guardrails(self):
        payload = self.payload
        lemma = payload["finite_hebbian_spectral_lemma"]
        self.assertIn("circulant", lemma["statement"])
        self.assertIn("not a strict nadsoliton theorem", lemma["proof_status"])
        impact = payload["impact_on_recent_highest_resonance_work"]
        self.assertIn("fifth-step/d5 orbit", impact["positive_impact"])
        self.assertIn("does not create the fifth-mode premise", impact["negative_result"])
        self.assertIn("no QW-2191 discharge", impact["selector_status"])
        ontology = payload["ontology_guardrail"]
        self.assertIn("primordial information", ontology["allowed_reading"])
        self.assertIn("not introduced as a separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("not as a sub-nadsoliton ontology", hard_limits)
        self.assertIn("Unbiased Hebbian training", hard_limits)
        self.assertIn("contiguous selects k=1", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
