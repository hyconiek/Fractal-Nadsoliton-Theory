import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3028_s1978_nadsoliton_information_to_classical_transition_foundation_lattice.py"
OUT = ROOT / "generated" / "p3028_s1978_nadsoliton_information_to_classical_transition_foundation_lattice.json"
MD = ROOT / "generated" / "p3028_s1978_nadsoliton_information_to_classical_transition_foundation_lattice.md"

class P3028InformationToClassicalFoundationLatticeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P3028_NADSOLITON_INFORMATION_TO_CLASSICAL_TRANSITION_FOUNDATION_LATTICE_NO_CLOSURE")
        self.assertIsNotNone(self.payload["input_hashes"]["P3013"])
        self.assertIsNotNone(self.payload["input_hashes"]["P3027"])

    def test_lattice_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["foundation_atom_count"], 5)
        self.assertEqual(cert["classical_readout_row_count"], 5)
        self.assertEqual(cert["accepted_readout_rows"], 0)
        self.assertEqual(cert["closure_profile_count"], 32)
        self.assertEqual(cert["accepting_profile_count"], 1)
        self.assertEqual(cert["current_closed_foundation_atoms"], 0)

    def test_rows_and_exports(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["object"], "NadsolitonInformationToClassicalTransition_FoundationObligationLattice")
        rows = obj["classical_readout_rows"]
        self.assertEqual([row["classical_row"] for row in rows], ["spacetime_geometry", "time", "matter_fields", "energy_hamiltonian", "observer_readout"])
        self.assertTrue(all(not row["accepted_as_classical_export"] for row in rows))
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3028/S1978", MD.read_text(encoding="utf-8"))
        self.assertIn("P3028/S1978", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3028/S1978", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3028", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
