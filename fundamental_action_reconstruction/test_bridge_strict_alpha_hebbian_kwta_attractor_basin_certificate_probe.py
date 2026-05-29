import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_kwta_attractor_basin_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_kwta_attractor_basin_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_kwta_attractor_basin_certificate_report.md"


class StrictAlphaHebbianKwtaAttractorBasinCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_finite_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_HEBBIAN_KWTA_ATTRACTOR_BASIN_CERTIFICATE_PROBE__NOT_A_THEOREM",
        )
        self.assertEqual(
            payload["status"],
            "conditional-d5-teacher-hebbian-dynamics-certificate-with-exact-tie-nonclosure",
        )
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["active_count"], 5)
        self.assertEqual(model["support_count"], 792)
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")

    def test_exact_weight_certificate(self):
        cert = self.payload["exact_weight_certificate"]
        self.assertEqual(
            cert["first_row_by_cyclic_distance"],
            ["0", "-25/12", "11/12", "-1/12", "-13/12", "23/12", "-25/12", "23/12", "-13/12", "-1/12", "11/12", "-25/12"],
        )
        self.assertTrue(cert["circulant_check"])
        self.assertTrue(cert["diagonal_zero_check"])
        self.assertIn("23/12", cert["positive_d5_neighbor_weight"])
        self.assertIn("-25/12", cert["negative_opposite_parity_weight"])

    def test_deterministic_kwta_reaches_d5_mostly_but_not_all(self):
        cert = self.payload["deterministic_lexicographic_kwta_certificate"]
        self.assertTrue(cert["all_teacher_patterns_fixed"])
        self.assertFalse(cert["all_initial_states_reach_d5"])
        self.assertEqual(cert["support_count"], 792)
        self.assertEqual(cert["reached_d5_count"], 786)
        self.assertEqual(cert["non_d5_attractor_count"], 6)
        self.assertEqual(cert["max_steps_to_d5"], 4)
        self.assertEqual(cert["steps_to_d5_histogram"], {"0": 12, "1": 444, "2": 289, "3": 33, "4": 8})
        self.assertEqual(len(cert["d5_basin_counts"]), 12)
        self.assertEqual(cert["non_d5_attractor_examples"][0]["initial"], [0, 1, 2, 3, 4])

    def test_set_valued_tie_branch_audit(self):
        audit = self.payload["set_valued_tie_branch_audit"]
        self.assertEqual(audit["support_count"], 792)
        self.assertEqual(audit["branch_count_histogram"], {"1": 624, "2": 120, "3": 24, "10": 12, "21": 12})
        self.assertEqual(audit["boundary_tie_support_count"], 168)
        self.assertEqual(audit["boundary_tie_size_histogram"], {"2": 120, "3": 24, "5": 12, "7": 12})
        self.assertEqual(audit["guaranteed_all_tie_branches_reach_d5_count"], 768)
        self.assertEqual(audit["guaranteed_all_tie_branches_reach_d5_fraction"], "768/792")
        self.assertEqual(audit["existential_some_tie_branch_reaches_d5_count"], 792)
        self.assertEqual(audit["existential_some_tie_branch_reaches_d5_fraction"], "792/792")
        self.assertEqual(audit["tie_sensitive_blocker_count"], 24)
        self.assertEqual(audit["tie_sensitive_blocker_fraction"], "24/792")

    def test_interpretation_and_guardrails(self):
        payload = self.payload
        impact = payload["impact_on_hebbian_hypothesis"]
        self.assertIn("786/792 states", impact["positive_gain"])
        self.assertIn("6/792 states enter non-d5", impact["honest_limit"])
        self.assertIn("only 768/792 states are guaranteed", impact["honest_limit"])
        self.assertIn("d5 teacher trace and tie selector", impact["meaning_for_highest_resonance"])
        self.assertIn("no QW-2191 discharge", impact["selector_status"])

        ontology = payload["ontology_guardrail"]
        self.assertIn("primordial information", ontology["allowed_reading"])
        self.assertIn("not a separate informational substrate", ontology["forbidden_reading"])

        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives the d5 teacher/self-record trace", hard_limits)
        self.assertIn("6/792 states enter non-d5", hard_limits)
        self.assertIn("24/792 states are tie-sensitive blockers", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
