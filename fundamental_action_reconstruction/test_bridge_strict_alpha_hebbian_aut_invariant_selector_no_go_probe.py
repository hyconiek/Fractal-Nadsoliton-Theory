import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_aut_invariant_selector_no_go_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_aut_invariant_selector_no_go_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_aut_invariant_selector_no_go_report.md"


class StrictAlphaHebbianAutInvariantSelectorNoGoProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_HEBBIAN_AUT_INVARIANT_SELECTOR_NO_GO_PROBE__NOT_A_THEOREM",
        )
        self.assertEqual(payload["status"], "finite-aut-invariant-selector-no-go-for-d5-unit-pair")
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["active_count"], 5)
        self.assertEqual(model["support_count"], 792)
        self.assertEqual(model["automorphism_units"], [1, 5, 7, 11])
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")

    def test_action_table_and_orbits(self):
        payload = self.payload
        action = payload["aut_action_table"]
        self.assertEqual(action["5"]["contiguous_step_1_11"], "fifth_step_d5_step_5_7")
        self.assertEqual(action["5"]["fifth_step_d5_step_5_7"], "contiguous_step_1_11")
        self.assertEqual(action["7"]["contiguous_step_1_11"], "fifth_step_d5_step_5_7")
        self.assertEqual(action["11"]["fifth_step_d5_step_5_7"], "fifth_step_d5_step_5_7")
        orbits = payload["orbit_partition"]
        self.assertEqual(orbits["orbit_count"], 2)
        self.assertEqual(orbits["unit_pair_orbit"], ["contiguous_step_1_11", "fifth_step_d5_step_5_7"])
        self.assertEqual(orbits["parity_orbit"], ["parity_minus_one_step_2_10"])

    def test_invariant_subset_enumeration(self):
        enum = self.payload["invariant_subset_enumeration"]
        self.assertEqual(enum["candidate_subset_count"], 8)
        self.assertEqual(enum["invariant_subset_count"], 4)
        self.assertEqual(
            enum["invariant_subsets"],
            [
                [],
                ["parity_minus_one_step_2_10"],
                ["contiguous_step_1_11", "fifth_step_d5_step_5_7"],
                ["contiguous_step_1_11", "parity_minus_one_step_2_10", "fifth_step_d5_step_5_7"],
            ],
        )

    def test_d5_selector_no_go_and_previous_selector(self):
        payload = self.payload
        no_go = payload["d5_selector_no_go_certificate"]
        self.assertEqual(no_go["d5_singleton"], ["fifth_step_d5_step_5_7"])
        self.assertFalse(no_go["d5_singleton_is_aut_invariant"])
        self.assertTrue(no_go["unit_pair_is_aut_invariant"])
        self.assertEqual(no_go["aut_invariant_singletons"], [["parity_minus_one_step_2_10"]])
        self.assertFalse(no_go["can_aut_invariant_selector_pick_d5_from_unit_pair"])
        self.assertIn("smallest Aut-invariant subset", no_go["reason"])

        previous = payload["previous_selector_classification"]
        self.assertEqual(previous["selected_subset"], ["fifth_step_d5_step_5_7"])
        self.assertFalse(previous["is_aut_invariant"])
        self.assertEqual(previous["selected_subset_images_by_unit"]["5"], ["contiguous_step_1_11"])
        self.assertIn("symmetry-breaking", previous["classification"])

    def test_proof_reading_and_guardrails(self):
        payload = self.payload
        proof = payload["proof_reading"]
        self.assertIn("union of the orbits", proof["finite_theorem"])
        self.assertIn("No Aut-invariant selector", proof["negative_result"])
        self.assertIn("Aut-breaking", proof["relation_to_previous_probe"])
        self.assertIn("orientation", proof["remaining_gap"])

        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives the d5 teacher/self-record trace", hard_limits)
        self.assertIn("No Aut(Z_12)-invariant selector", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
