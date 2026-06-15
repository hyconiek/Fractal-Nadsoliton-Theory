import json
import subprocess
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "fundamental_action_reconstruction"
SCRIPT = ROOT / "p2771_s1721_geometric_self_coupling_operator_witness.py"
OUT = ROOT / "generated" / "p2771_s1721_geometric_self_coupling_operator_witness.json"
MD = ROOT / "generated" / "p2771_s1721_geometric_self_coupling_operator_witness.md"


class P2771GeometricSelfCouplingOperatorWitnessTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run(["python", str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2771_FINITE_GEOMETRIC_SELF_COUPLING_OPERATOR_WITNESS_BOUNDED_NO_GO_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2770"], "P2770_KERNEL_CHARACTERISTIC_EXPRESSIVITY_AUDIT_NO_FULL_EXPRESSION_NO_CLOSURE")
        self.assertEqual(self.payload["audited_object"], "selector-free finite cyclic radial geometric self-coupling operator C_geo_N[K]")

    def test_c_geo_witness_fails_for_both_kernels(self):
        witness = self.payload["geometric_self_coupling_witness"]
        self.assertIn("C_geo_N[K]", witness["operator"])
        self.assertEqual(witness["row_count"], 2)
        self.assertFalse(witness["all_kernels_pass_scalar_eigenclosure"])
        self.assertEqual(set(witness["failed_kernels"]), {"K_legacy_ont", "K_strict_gate"})
        self.assertGreater(witness["max_relative_l2_residual"], 1e-9)
        for row in witness["rows"]:
            self.assertEqual(row["shell_order"], [0, 1, 2, 3, 4, 5, 6])
            self.assertFalse(row["passes_scalar_eigenclosure_tolerance_1e_minus_9"])
            self.assertGreater(row["relative_l2_residual"], 1e-9)
            self.assertGreater(row["max_abs_residual"], 1e-9)

    def test_acceptance_blocks_geometric_closure_and_ltotal(self):
        acceptance = self.payload["acceptance_matrix"]
        self.assertTrue(acceptance["facts"]["explicit_c_geo_candidate_supplied"])
        self.assertFalse(acceptance["facts"]["all_kernel_samples_scalar_eigenclosed"])
        self.assertFalse(acceptance["facts"]["candidate_is_unique_or_ontologically_sourced"])
        self.assertFalse(acceptance["accepted_as_geometric_self_coupling_theorem"])
        self.assertFalse(acceptance["accepted_as_kernel_full_expression_theorem"])
        self.assertFalse(acceptance["accepted_as_learning_update_law"])
        self.assertFalse(acceptance["accepted_as_ltotal_or_toe_promotion"])
        self.assertIn("all_kernel_samples_scalar_eigenclosed", acceptance["missing_criteria"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_documentation_and_recommendation(self):
        self.assertIn("self-learning kernel-parameter update law", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2697-P2771", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2771/S1721", MD.read_text(encoding="utf-8"))
        self.assertIn("P2771/S1721", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2771/S1721", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2771", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
