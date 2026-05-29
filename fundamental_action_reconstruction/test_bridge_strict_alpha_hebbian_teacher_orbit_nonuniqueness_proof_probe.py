import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_teacher_orbit_nonuniqueness_proof_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_teacher_orbit_nonuniqueness_proof_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_teacher_orbit_nonuniqueness_proof_report.md"


class StrictAlphaHebbianTeacherOrbitNonuniquenessProofProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_HEBBIAN_TEACHER_ORBIT_NONUNIQUENESS_PROOF_PROBE__NOT_A_THEOREM",
        )
        self.assertEqual(payload["status"], "finite-hebbian-associative-memory-nonuniqueness-not-d5-source")
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["active_count"], 5)
        self.assertEqual(model["support_count"], 792)
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")

    def test_each_teacher_class_self_maximizes(self):
        certificates = self.payload["orbit_class_certificates"]
        self.assertEqual(set(certificates), {
            "contiguous_step_1_11",
            "parity_minus_one_step_2_10",
            "fifth_step_d5_step_5_7",
        })
        for certificate in certificates.values():
            self.assertTrue(certificate["step_orbits_equal"])
            self.assertEqual(certificate["teacher_support_count"], 12)
            self.assertEqual(certificate["global_maximizer_count"], 12)
            self.assertTrue(certificate["global_maximizers_equal_teacher_orbit"])
            self.assertEqual(certificate["non_teacher_global_maximizer_count"], 0)
            self.assertEqual(certificate["histogram_class_count"], 35)
            self.assertTrue(certificate["direct_energy_equals_histogram_replay_for_all_supports"])
            self.assertEqual(certificate["replay_failure_count"], 0)

    def test_class_specific_certificates(self):
        certificates = self.payload["orbit_class_certificates"]
        contiguous = certificates["contiguous_step_1_11"]
        parity = certificates["parity_minus_one_step_2_10"]
        d5 = certificates["fifth_step_d5_step_5_7"]

        self.assertEqual(contiguous["maximum_energy"], "55/3")
        self.assertEqual(contiguous["top_distance_histogram"]["distance_histogram_d1_to_d6"], [4, 3, 2, 1, 0, 0])
        self.assertEqual(parity["maximum_energy"], "115/3")
        self.assertEqual(parity["top_distance_histogram"]["distance_histogram_d1_to_d6"], [0, 4, 0, 4, 0, 2])
        self.assertEqual(d5["maximum_energy"], "55/3")
        self.assertEqual(d5["top_distance_histogram"]["distance_histogram_d1_to_d6"], [0, 3, 2, 1, 4, 0])
        self.assertEqual(d5["first_row_by_cyclic_distance"], ["0", "-25/12", "11/12", "-1/12", "-13/12", "23/12", "-25/12", "23/12", "-13/12", "-1/12", "11/12", "-25/12"])

    def test_cross_orbit_nonuniqueness(self):
        cross = self.payload["cross_orbit_nonuniqueness_certificate"]
        self.assertTrue(cross["all_classes_self_maximize"])
        self.assertEqual(
            cross["non_d5_self_maximizing_classes"],
            ["contiguous_step_1_11", "parity_minus_one_step_2_10"],
        )
        self.assertEqual(cross["d5_maximum_energy"], "55/3")
        self.assertEqual(cross["d5_global_maximizer_count"], 12)
        self.assertIn("does not decide which teacher orbit", cross["honest_obstruction"])

    def test_proof_reading_and_guardrails(self):
        payload = self.payload
        proof = payload["proof_reading"]
        self.assertIn("each admissible", proof["finite_theorem"])
        self.assertIn("not source-selective", proof["negative_consequence_for_d5_source"])
        self.assertIn("associative-memory certificate", proof["relation_to_previous_probe"])
        self.assertIn("strict-side source", proof["remaining_gap"])

        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives the d5 teacher/self-record trace", hard_limits)
        self.assertIn("nonunique across teacher orbits", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
