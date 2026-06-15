import json
import subprocess
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "fundamental_action_reconstruction"
SCRIPT = ROOT / "p2772_s1722_self_learning_kernel_update_law_stationarity_witness.py"
OUT = ROOT / "generated" / "p2772_s1722_self_learning_kernel_update_law_stationarity_witness.json"
MD = ROOT / "generated" / "p2772_s1722_self_learning_kernel_update_law_stationarity_witness.md"


class P2772SelfLearningKernelUpdateLawStationarityWitnessTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run(["python", str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2772_SELF_LEARNING_KERNEL_UPDATE_LAW_STATIONARITY_WITNESS_BOUNDED_NO_GO_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2771"], "P2771_FINITE_GEOMETRIC_SELF_COUPLING_OPERATOR_WITNESS_BOUNDED_NO_GO_NO_CLOSURE")
        self.assertEqual(self.payload["audited_object"], "finite gradient self-learning update law for the P2771 C_geo residual loss")

    def test_update_witness_has_nonzero_gradients_for_both_kernels(self):
        witness = self.payload["self_learning_update_witness"]
        self.assertEqual(witness["candidate_update_law"], "theta_{t+1}=theta_t - lr * grad_theta L_geo(theta_t)")
        self.assertEqual(witness["row_count"], 2)
        self.assertFalse(witness["all_current_tuples_stationary"])
        self.assertEqual(set(witness["failed_stationary_kernels"]), {"K_legacy_ont", "K_strict_gate"})
        self.assertGreater(witness["min_gradient_norm"], 1e-9)
        for row in witness["rows"]:
            self.assertFalse(row["passes_stationary_tolerance"])
            self.assertGreater(row["gradient_norm"], 1e-9)
            self.assertGreater(len(row["nonzero_update_parameters"]), 0)
            self.assertIn("L_geo", row["loss_definition"])

    def test_acceptance_blocks_learning_closure_and_ltotal(self):
        acceptance = self.payload["acceptance_matrix"]
        self.assertTrue(acceptance["facts"]["explicit_update_law_supplied"])
        self.assertFalse(acceptance["facts"]["all_current_kernel_tuples_stationary"])
        self.assertFalse(acceptance["facts"]["candidate_loss_ontologically_sourced"])
        self.assertFalse(acceptance["accepted_as_self_learning_kernel_theorem"])
        self.assertFalse(acceptance["accepted_as_learned_stationary_fixed_point"])
        self.assertFalse(acceptance["accepted_as_geometric_self_coupling_theorem"])
        self.assertFalse(acceptance["accepted_as_ltotal_or_toe_promotion"])
        self.assertIn("all_current_kernel_tuples_stationary", acceptance["missing_criteria"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_documentation_and_recommendation(self):
        self.assertIn("ontologically sourced learning functional", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2697-P2772", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2772/S1722", MD.read_text(encoding="utf-8"))
        self.assertIn("P2772/S1722", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2772/S1722", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2772", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
