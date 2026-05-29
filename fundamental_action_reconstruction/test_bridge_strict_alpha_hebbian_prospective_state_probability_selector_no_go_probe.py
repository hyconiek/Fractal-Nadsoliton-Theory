import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_prospective_state_probability_selector_no_go_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_prospective_state_probability_selector_no_go_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_prospective_state_probability_selector_no_go_report.md"


class StrictAlphaHebbianProspectiveStateProbabilitySelectorNoGoProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_HEBBIAN_PROSPECTIVE_STATE_PROBABILITY_SELECTOR_NO_GO_PROBE__NOT_A_THEOREM",
        )
        self.assertEqual(payload["status"], "aut-invariant-prospective-state-probability-cannot-break-unit-mirror-tie")
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["support_count"], 792)
        self.assertEqual(model["histogram_class_count"], 35)
        self.assertEqual(model["automorphism_units"], [1, 5, 7, 11])
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")
        self.assertEqual([axis["name"] for axis in model["survivor_axes"]], ["A1_k1_contiguous", "A5_k5_d5"])

    def test_unit_action_and_future_predicates(self):
        cert = self.payload["unit_action_certificate"]
        self.assertEqual(cert["swapping_units"], [5, 7])
        self.assertEqual(cert["preserving_units"], [1, 11])
        by_unit = {row["unit"]: row for row in cert["action_table"]}
        self.assertEqual(by_unit[5]["maps_A1_to"], "A5_k5_d5")
        self.assertEqual(by_unit[5]["maps_A5_to"], "A1_k1_contiguous")

        predicates = self.payload["future_predicate_enumeration"]
        self.assertEqual(predicates["predicate_count"], 4)
        self.assertEqual(predicates["full_aut_invariant_predicate_names"], ["none", "both"])
        self.assertFalse(predicates["full_aut_singleton_d5_future_predicate_exists"])

    def test_prospective_prior_grid_audit(self):
        audit = self.payload["prospective_prior_grid_audit"]
        self.assertEqual(audit["denominator"], 12)
        self.assertEqual(audit["prior_count"], 13)
        self.assertEqual(audit["full_aut_invariant_prior_count"], 1)
        self.assertEqual(audit["d5_dominant_prior_count"], 6)
        self.assertEqual(audit["full_aut_invariant_d5_dominant_prior_count"], 0)
        invariant_rows = [row for row in audit["rows"] if row["full_aut_invariant"]]
        self.assertEqual(invariant_rows, [{"P_A1": "1/2", "P_A5": "1/2", "full_aut_invariant": True, "d5_dominant": False, "imports_unit_orientation_bit": False}])

    def test_prospective_transition_grid_audit(self):
        audit = self.payload["prospective_transition_grid_audit"]
        self.assertEqual(audit["denominator"], 12)
        self.assertEqual(audit["kernel_count"], 169)
        self.assertEqual(audit["full_aut_equivariant_kernel_count"], 13)
        self.assertEqual(audit["d5_dominant_from_uniform_prior_count"], 78)
        self.assertEqual(audit["full_aut_equivariant_d5_dominant_from_uniform_prior_count"], 0)
        self.assertIn("[[a, 1-a], [1-a, a]]", audit["full_aut_equivariant_family"])
        self.assertTrue(all(not row["d5_dominant_from_uniform_prior"] for row in audit["sample_equivariant_kernels"]))
        self.assertTrue(all(row["imports_unit_orientation_bit_if_d5_dominant"] for row in audit["sample_non_equivariant_d5_dominant_kernels"]))

    def test_deterministic_maps_interpretation_and_guardrails(self):
        maps = self.payload["deterministic_future_map_enumeration"]
        self.assertEqual(maps["map_count"], 4)
        self.assertEqual(maps["full_aut_equivariant_map_count"], 2)
        self.assertTrue(maps["equivariant_maps_do_not_create_new_orientation_from_unoriented_orbit"])

        interpretation = self.payload["selector_interpretation"]
        self.assertIn("prospective self-record", interpretation["faith_like_hypothesis_tested"])
        self.assertIn("equal probability", interpretation["negative_result"])
        self.assertIn("one-bit unit-axis", interpretation["conditional_positive_result"])
        self.assertIn("highest-future-resonance", interpretation["highest_resonance_reading"])

        ontology = self.payload["ontology_guardrail"]
        self.assertIn("nadsoliton itself is informational", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives a faith-like prospective self-record", hard_limits)
        self.assertIn("Full Aut(Z_12)-invariant prospective probabilities", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
