import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_energy_refined_kwta_selector_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_energy_refined_kwta_selector_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_energy_refined_kwta_selector_report.md"


class StrictAlphaHebbianEnergyRefinedKwtaSelectorProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_HEBBIAN_ENERGY_REFINED_KWTA_SELECTOR_PROBE__NOT_A_THEOREM",
        )
        self.assertEqual(
            payload["status"],
            "conditional-energy-refined-hebbian-selector-closes-finite-d5-basin-not-strict-source",
        )
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["active_count"], 5)
        self.assertEqual(model["support_count"], 792)
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")
        self.assertIn("maximizing learned energy", model["tie_refinement"])

    def test_exact_weight_replay_and_branch_refinement(self):
        payload = self.payload
        weight = payload["exact_weight_replay"]
        self.assertEqual(
            weight["first_row_by_cyclic_distance"],
            ["0", "-25/12", "11/12", "-1/12", "-13/12", "23/12", "-25/12", "23/12", "-13/12", "-1/12", "11/12", "-25/12"],
        )
        self.assertTrue(weight["circulant_check"])
        self.assertTrue(weight["diagonal_zero_check"])

        audit = payload["branch_refinement_audit"]
        self.assertEqual(audit["raw_kwta_branch_count_histogram"], {"1": 624, "2": 120, "3": 24, "10": 12, "21": 12})
        self.assertEqual(audit["energy_refined_branch_count_histogram"], {"1": 720, "2": 60, "4": 12})
        self.assertEqual(audit["supports_changed_by_energy_refinement_count"], 120)
        self.assertEqual(audit["supports_still_multibranch_after_energy_refinement_count"], 72)

    def test_all_retained_branches_close_to_d5(self):
        cert = self.payload["energy_refined_all_branch_closure_certificate"]
        self.assertEqual(cert["guaranteed_all_retained_branches_reach_d5_count"], 792)
        self.assertEqual(cert["guaranteed_all_retained_branches_reach_d5_fraction"], "792/792")
        self.assertEqual(cert["closure_layer_histogram"], {"0": 12, "1": 615, "2": 156, "3": 9})
        self.assertEqual(cert["max_closure_layer"], 3)
        self.assertEqual(cert["unclosed_count"], 0)

    def test_deterministic_replay_and_energy_ascent(self):
        payload = self.payload
        replay = payload["deterministic_replay_with_lex_only_inside_equal_energy_ties"]
        self.assertTrue(replay["all_initial_states_reach_d5"])
        self.assertEqual(replay["reached_d5_count"], 792)
        self.assertEqual(replay["failure_count"], 0)
        self.assertEqual(replay["steps_to_d5_histogram"], {"0": 12, "1": 480, "2": 264, "3": 36})
        self.assertEqual(len(replay["d5_basin_counts"]), 12)

        energy = payload["energy_ascent_certificate"]
        self.assertTrue(energy["strict_positive_ascent_for_every_retained_non_d5_branch"])
        self.assertEqual(energy["minimum_non_d5_delta"], "4")
        self.assertEqual(energy["violation_count"], 0)
        self.assertEqual(energy["delta_histogram_over_retained_branches"]["0"], 12)
        self.assertEqual(energy["delta_histogram_over_retained_branches"]["28"], 132)

    def test_interpretation_and_guardrails(self):
        payload = self.payload
        impact = payload["impact_on_previous_kwta_nonclosure"]
        self.assertIn("6/792 deterministic non-d5 cycles", impact["positive_gain"])
        self.assertIn("24/792 raw tie-sensitive blockers", impact["positive_gain"])
        self.assertIn("neither premise is derived", impact["remaining_limit"])
        self.assertIn("not a QW-2191 discharge", impact["selector_status"])

        ontology = payload["ontology_guardrail"]
        self.assertIn("finite self-recorded resonance patterns", ontology["allowed_reading"])
        self.assertIn("not a separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives the d5 teacher/self-record trace", hard_limits)
        self.assertIn("No theorem derives the centered zero-self Hebbian rule", hard_limits)
        self.assertIn("not strict-core selector closure", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
