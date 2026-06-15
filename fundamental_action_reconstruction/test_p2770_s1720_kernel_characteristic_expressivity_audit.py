import json
import subprocess
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "fundamental_action_reconstruction"
SCRIPT = ROOT / "p2770_s1720_kernel_characteristic_expressivity_audit.py"
OUT = ROOT / "generated" / "p2770_s1720_kernel_characteristic_expressivity_audit.json"
MD = ROOT / "generated" / "p2770_s1720_kernel_characteristic_expressivity_audit.md"


class P2770KernelCharacteristicExpressivityAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run(["python", str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2770_KERNEL_CHARACTERISTIC_EXPRESSIVITY_AUDIT_NO_FULL_EXPRESSION_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2769"], "P2769_COMBINED_NORMALIZATION_ORBIT_TRANSITIVITY_NO_GO_NO_CLOSURE")
        self.assertIn("K_legacy_ont", self.payload["kernels"])
        self.assertIn("K_strict_gate", self.payload["kernels"])

    def test_formula_matrix_has_bounded_solitonic_coverage_only(self):
        matrix = self.payload["formula_expressivity_matrix"]
        self.assertEqual(matrix["row_count"], 5)
        self.assertEqual(matrix["kernel_count"], 2)
        self.assertFalse(matrix["all_characteristics_formula_level_pass"])
        self.assertIn("self_learning_neural_network_character", matrix["failed_characteristics"])
        self.assertIn("geometric_self_coupling", matrix["failed_characteristics"])
        solitonic = next(row for row in matrix["rows"] if row["characteristic_id"] == "solitonic_state")
        self.assertTrue(solitonic["any_kernel_formula_level_pass"])
        self.assertTrue(solitonic["kernels"]["K_legacy_ont"]["formula_level_pass"])
        self.assertTrue(solitonic["kernels"]["K_strict_gate"]["formula_level_pass"])

    def test_finite_feature_rank_witness_records_defects(self):
        witness = self.payload["finite_feature_rank_witness"]
        self.assertEqual(witness["combined_formula_coverage_count"], 1)
        self.assertEqual(witness["required_coverage_count"], 5)
        self.assertEqual(witness["coverage_defect_count"], 4)
        self.assertIn("single_primordial_foundation", witness["coverage_defect_ids"])
        self.assertIn("complete_kernel_expression", witness["coverage_defect_ids"])
        self.assertIn("solitonic", witness["finite_conclusion"])

    def test_acceptance_blocks_full_expression_and_closure(self):
        acceptance = self.payload["acceptance_matrix"]
        self.assertFalse(acceptance["accepted_as_full_kernel_characterization"])
        self.assertFalse(acceptance["accepted_as_self_learning_kernel_theorem"])
        self.assertFalse(acceptance["accepted_as_geometric_self_coupling_kernel_theorem"])
        self.assertFalse(acceptance["accepted_as_legacy_strict_bridge_closure"])
        self.assertFalse(acceptance["accepted_as_ltotal_or_toe_promotion"])
        self.assertIn("formula_covers_all_required_characteristics", acceptance["missing_criteria"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_documentation_and_recommendation(self):
        self.assertIn("C_geo[K]", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2697-P2770", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2770/S1720", MD.read_text(encoding="utf-8"))
        self.assertIn("P2770/S1720", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2770/S1720", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2770", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
