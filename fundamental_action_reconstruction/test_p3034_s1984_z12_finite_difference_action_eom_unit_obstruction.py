import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3034_s1984_z12_finite_difference_action_eom_unit_obstruction.py"
OUT = ROOT / "generated" / "p3034_s1984_z12_finite_difference_action_eom_unit_obstruction.json"
MD = ROOT / "generated" / "p3034_s1984_z12_finite_difference_action_eom_unit_obstruction.md"

class P3034Z12FiniteDifferenceActionEOMUnitTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P3034_Z12_FINITE_DIFFERENCE_ACTION_EOM_UNIT_OBSTRUCTION_NO_LTOTAL_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P3028"])
        self.assertIsNotNone(self.payload["input_hashes"]["P3033"])

    def test_finite_certificate(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["receiver_rows"], 3)
        self.assertEqual(cert["unit_bearing_receiver_rows"], 0)
        self.assertEqual(cert["action_rescaling_exponent"], 2.0)
        self.assertEqual(cert["residual_rescaling_exponent"], 1.0)
        self.assertFalse(cert["unit_bearing_action_eom_hamiltonian_exported"])

    def test_obligations_and_exports(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["object"], "Z12FiniteDifferenceActionEOMHamiltonian_UnitObstructionMatrix")
        obligations = {row["obligation"]: row["satisfied"] for row in obj["proof_obligations"]}
        self.assertTrue(obligations["attacks_single_P3028_foundation_atom"])
        self.assertTrue(obligations["explicit_formal_action_receiver"])
        self.assertTrue(obligations["finite_eom_residual_computable"])
        self.assertTrue(obligations["cyclic_shift_covariant_receiver"])
        self.assertFalse(obligations["physical_action_unit_source"])
        self.assertFalse(obligations["field_provenance_and_boundary_map"])
        self.assertFalse(obligations["physical_clock_or_hamiltonian_unit"])
        self.assertFalse(obligations["nonproxy_continuum_eom_lift"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3034/S1984", MD.read_text(encoding="utf-8"))
        self.assertIn("P3034/S1984", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3034/S1984", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3034", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
