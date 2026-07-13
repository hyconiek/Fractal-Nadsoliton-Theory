import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3142_s2092_axiom_selector_field_variational_lift.py"
OUT = ROOT / "generated" / "p3142_s2092_axiom_selector_field_variational_lift.json"
MD = ROOT / "generated" / "p3142_s2092_axiom_selector_field_variational_lift.md"


class P3142AxiomSelectorFieldVariationalLiftTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_inputs_and_counts(self):
        self.assertEqual(self.payload["status"], "P3142_AXIOM_SELECTOR_FIELD_VARIATIONAL_LIFT_BOUNDED_NON_STRICT")
        self.assertTrue(all(self.payload["input_hashes"].values()))
        counts = self.payload["finite_theorem"]["finite_counts"]
        self.assertEqual(counts["parameter_triples_tested"], 27)
        self.assertEqual(counts["stationary_rows"], 27)
        self.assertEqual(counts["positive_hessian_rows"], 27)
        self.assertEqual(counts["strict_weight_source_rows"], 0)
        self.assertEqual(counts["unit_bearing_measure_rows"], 0)

    def test_symbolic_variation(self):
        variation = self.payload["symbolic_variation"]
        self.assertEqual(variation["gradient_at_selected"], ["0", "0"])
        self.assertTrue(variation["positive_definite_for_positive_parameters"])
        self.assertIn("theta", variation["variables"])
        self.assertIn("s", variation["variables"])

    def test_parameter_rows_and_gates(self):
        self.assertEqual(len(self.payload["parameter_rows"]), 27)
        self.assertTrue(all(row["gradient_at_selected"] == ["0", "0"] for row in self.payload["parameter_rows"]))
        self.assertTrue(all(min(row["positive_eigenvalues"]) > 0 for row in self.payload["parameter_rows"]))
        gates = {row["gate"]: row["passed"] for row in self.payload["gate_rows"]}
        self.assertTrue(gates["local_variational_derivative"])
        self.assertTrue(gates["selected_point_stationary"])
        self.assertTrue(gates["local_positive_hessian"])
        self.assertFalse(gates["strict_source_weights"])
        self.assertFalse(gates["global_Z12_field_chart"])
        self.assertFalse(gates["unit_bearing_L_total_ToE"])

    def test_docs_updated(self):
        self.assertIn("P3142/S2092", MD.read_text(encoding="utf-8"))
        self.assertIn("P3142/S2092", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3142/S2092", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3142", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
