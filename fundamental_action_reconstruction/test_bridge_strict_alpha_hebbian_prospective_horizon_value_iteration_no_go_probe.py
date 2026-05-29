import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_prospective_horizon_value_iteration_no_go_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_prospective_horizon_value_iteration_no_go_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_prospective_horizon_value_iteration_no_go_report.md"


class StrictAlphaHebbianProspectiveHorizonValueIterationNoGoProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_finite_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_HEBBIAN_PROSPECTIVE_HORIZON_VALUE_ITERATION_NO_GO_PROBE__NOT_A_THEOREM",
        )
        self.assertEqual(
            payload["status"],
            "finite-horizon-aut-equivariant-lookahead-cannot-create-d5-selector-from-invariant-future-value",
        )
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["support_count"], 792)
        self.assertEqual(model["histogram_class_count"], 35)
        self.assertEqual(model["automorphism_units"], [1, 5, 7, 11])
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")
        self.assertEqual([axis["name"] for axis in model["survivor_axes"]], ["A1_k1_contiguous", "A5_k5_d5"])

    def test_equivariant_transition_family_and_reward_grid(self):
        family = self.payload["equivariant_transition_family"]
        self.assertEqual(family["denominator"], 12)
        self.assertEqual(family["kernel_count"], 13)
        self.assertEqual(family["closed_form"], "T_a = [[a, 1-a], [1-a, a]]")
        self.assertTrue(family["commutes_with_unit_mirror_swap"])

        audit = self.payload["future_reward_grid_audit"]
        self.assertEqual(audit["reward_count"], 169)
        self.assertEqual(audit["full_aut_invariant_reward_count"], 13)
        self.assertEqual(audit["d5_dominant_reward_count"], 78)
        self.assertEqual(audit["full_aut_invariant_d5_dominant_reward_count"], 0)
        self.assertTrue(all(row["tied"] for row in audit["sample_invariant_rewards"]))
        self.assertTrue(all(row["d5_dominant"] for row in audit["sample_d5_dominant_rewards"]))

    def test_fixed_kernel_finite_horizon_no_go(self):
        audit = self.payload["fixed_kernel_finite_horizon_audit"]
        self.assertEqual(audit["max_horizon"], 8)
        self.assertTrue(audit["all_invariant_rewards_remain_tied"])
        self.assertEqual(audit["invariant_reward_value_tie_violations"], 0)
        self.assertEqual([row["horizon"] for row in audit["rows"]], list(range(9)))
        for row in audit["rows"]:
            self.assertEqual(row["cases_checked"], 169)
            self.assertEqual(row["tied_value_count"], 169)
            self.assertEqual(row["d5_dominant_value_count"], 0)

    def test_controlled_possible_future_and_proof_certificate(self):
        controlled = self.payload["controlled_possible_future_audit"]
        self.assertTrue(controlled["all_controlled_invariant_values_remain_tied"])
        self.assertTrue(all(row["tied"] for row in controlled["invariant_terminal_bellman_trace"]))
        labelled_samples = controlled["d5_labelled_terminal_samples"]
        self.assertIn("d5 terminal value is non-invariant", labelled_samples[0]["interpretation"])
        self.assertIn("does not select a unique current axis", labelled_samples[-1]["interpretation"])
        self.assertTrue(labelled_samples[-1]["bellman_max_values_by_horizon"][1]["tied"])

        proof = self.payload["exact_proof_certificate"]
        self.assertEqual(proof["equivariance_condition"], "T J = J T")
        self.assertIn("T_a", proof["two_state_solution"])
        self.assertIn("finite horizon", proof["finite_horizon_induction"])

    def test_interpretation_and_guardrails(self):
        interpretation = self.payload["selector_interpretation"]
        self.assertIn("highest possible future resonance", interpretation["question_tested"])
        self.assertIn("every finite horizon remains tied", interpretation["negative_result"])
        self.assertIn("one-bit unit-axis", interpretation["conditional_positive_result"])
        self.assertIn("not a theorem deriving", interpretation["honest_limit"])

        ontology = self.payload["ontology_guardrail"]
        self.assertIn("nadsoliton itself", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives a faith-like prospective value", hard_limits)
        self.assertIn("Full Aut(Z_12)-equivariant finite-horizon lookahead", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
