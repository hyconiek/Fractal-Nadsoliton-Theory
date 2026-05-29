import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_distance5_pair_count_core_selector_proof_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_distance5_pair_count_core_selector_proof_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_distance5_pair_count_core_selector_proof_report.md"


class StrictAlphaHebbianDistance5PairCountCoreSelectorProofProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_HEBBIAN_DISTANCE5_PAIR_COUNT_CORE_SELECTOR_PROOF_PROBE__NOT_A_THEOREM",
        )
        self.assertEqual(
            payload["status"],
            "distance5-pair-count-uniquely-selects-d5-conditionally-not-origin-theorem",
        )
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["active_count"], 5)
        self.assertEqual(model["support_count"], 792)
        self.assertEqual(model["target_distance"], 5)
        self.assertEqual(model["teacher_orbit_size"], 12)
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")

    def test_distance5_graph_certificate(self):
        graph = self.payload["distance5_graph_certificate"]
        self.assertEqual(graph["edge_count"], 12)
        self.assertEqual(graph["gcd_step_ring"], 1)
        self.assertTrue(graph["is_single_cycle"])
        self.assertEqual(graph["cycle_order_by_adding_5"], [0, 5, 10, 3, 8, 1, 6, 11, 4, 9, 2, 7])
        self.assertEqual(graph["path_support_count"], 12)
        self.assertTrue(graph["path_supports_equal_teacher_orbit"])
        self.assertIn("5-vertex induced subgraph of a cycle", graph["graph_bound"])

    def test_pair_count_selector_certificate(self):
        certificate = self.payload["pair_count_selector_certificate"]
        self.assertEqual(certificate["max_distance5_pair_count"], 4)
        self.assertEqual(certificate["maximizer_count"], 12)
        self.assertTrue(certificate["maximizers_equal_d5_teacher_orbit"])
        self.assertTrue(certificate["maximizers_equal_distance5_cycle_paths"])
        self.assertEqual(certificate["count_distribution"], {"0": 36, "1": 240, "2": 360, "3": 144, "4": 12})
        self.assertEqual(certificate["closest_nonmax_count"], 3)
        self.assertEqual(certificate["closest_nonmax_support_count"], 144)
        self.assertEqual(certificate["closest_nonmax_examples"][0], [0, 1, 2, 5, 7])
        self.assertEqual(certificate["objective_unit_stabilizer"], [1, 11])
        self.assertTrue(certificate["objective_stabilizer_equals_required_subgroup"])

    def test_relation_interpretation_and_guardrails(self):
        payload = self.payload
        relation = payload["relation_to_previous_score"]
        self.assertEqual(relation["previous_score_formula"], "C(h)=3*h_2 + 2*h_3 + h_4 + 4*h_5")
        self.assertEqual(relation["core_score_coefficients_d1_to_d6"], [0, 0, 0, 0, 1, 0])
        self.assertIn("h_5 term alone", relation["finite_gain"])

        interpretation = payload["candidate_source_interpretation"]
        self.assertIn("five vertices in the distance-5 cycle", interpretation["finite_gain"])
        self.assertIn("distance-5 pair-count maximization", interpretation["conditional_positive_result"])
        self.assertIn("does not derive why distance 5", interpretation["honest_limit"])

        ontology = payload["ontology_guardrail"]
        self.assertIn("nadsoliton itself", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives the d5 teacher/self-record trace", hard_limits)
        self.assertIn("No theorem derives a Hebbian learning law", hard_limits)
        self.assertIn("No theorem derives the distance-5 pair-count objective", hard_limits)
        self.assertIn("conditional on supplying distance 5", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
