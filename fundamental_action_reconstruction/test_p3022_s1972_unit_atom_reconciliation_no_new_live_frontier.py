import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3022_s1972_unit_atom_reconciliation_no_new_live_frontier.py"
OUT = ROOT / "generated" / "p3022_s1972_unit_atom_reconciliation_no_new_live_frontier.json"
MD = ROOT / "generated" / "p3022_s1972_unit_atom_reconciliation_no_new_live_frontier.md"

class P3022UnitAtomReconciliationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P3022_UNIT_ATOM_RECONCILIATION_NO_NEW_LIVE_FRONTIER_NO_CLOSURE")
        self.assertEqual(set(self.payload["input_hashes"].keys()), {"P3017", "P3018", "P3019", "P3020", "P3021"})
        for value in self.payload["input_hashes"].values():
            self.assertIsNotNone(value)

    def test_five_atom_lattice_no_current_acceptance(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["unit_atom_count"], 5)
        self.assertEqual(cert["closed_source_atom_count"], 0)
        self.assertEqual(cert["closure_profile_count"], 32)
        self.assertEqual(cert["accepting_profile_count"], 1)
        self.assertFalse(cert["current_profile_accepts"])
        self.assertEqual(cert["new_live_frontier_count"], 0)

    def test_rows_and_negative_exports(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["object"], "TimeObservableUnitAtomReconciliation_NoNewLiveFrontierCertificate")
        atoms = {row["atom"]: row for row in obj["unit_atom_rows"]}
        self.assertIn("formal_time_observable_action_eom", atoms)
        self.assertIn("lambda_action_unit_normalization", atoms)
        self.assertIn("observable_unit_readout", atoms)
        self.assertIn("clock_unit_theorem", atoms)
        self.assertIn("action_quantum_reference_cell_source", atoms)
        self.assertTrue(all(not row["source_closed"] for row in atoms.values()))
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3022/S1972", MD.read_text(encoding="utf-8"))
        self.assertIn("P3022/S1972", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3022/S1972", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3022", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
