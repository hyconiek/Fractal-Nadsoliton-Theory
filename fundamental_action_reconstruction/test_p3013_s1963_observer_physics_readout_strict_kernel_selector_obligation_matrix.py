import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3013_s1963_observer_physics_readout_strict_kernel_selector_obligation_matrix.py"
OUT = ROOT / "generated" / "p3013_s1963_observer_physics_readout_strict_kernel_selector_obligation_matrix.json"
MD = ROOT / "generated" / "p3013_s1963_observer_physics_readout_strict_kernel_selector_obligation_matrix.md"

class P3013ObserverPhysicsReadoutTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P3013_OBSERVER_PHYSICS_READOUT_STRICT_KERNEL_SELECTOR_OBLIGATION_MATRIX_NO_CLOSURE")
        self.assertIsNotNone(self.payload["input_hashes"]["P3011"])
        self.assertIsNotNone(self.payload["input_hashes"]["P3012"])

    def test_readout_certificate(self):
        cert = self.payload["readout_certificate"]
        self.assertEqual(cert["physics_rows"], ["spacetime_geometry", "time", "matter", "energy", "observer_readout"])
        self.assertEqual(cert["row_count"], 5)
        self.assertEqual(cert["accepted_row_count"], 0)
        self.assertEqual(cert["global_acceptance_profile_count"], 32)
        self.assertEqual(cert["global_accepting_profile_count"], 1)
        self.assertEqual(cert["blocked_atom_union"], ["eom_hamiltonian", "observable_generator", "selector_source", "unit_action"])

    def test_rows_and_nonpromotion(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["object"], "ObserverPhysicsReadout_StrictKernelSelector_ObligationMatrix")
        for row in obj["rows"]:
            self.assertTrue(row["current_atoms"]["kernel_profile"])
            self.assertFalse(row["accepted_as_observed_physics_export"])
            self.assertIn("selector_source", row["blocked_atoms"])
            self.assertIn("observable_generator", row["blocked_atoms"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3013/S1963", MD.read_text(encoding="utf-8"))
        self.assertIn("P3013/S1963", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3013/S1963", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3013", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
