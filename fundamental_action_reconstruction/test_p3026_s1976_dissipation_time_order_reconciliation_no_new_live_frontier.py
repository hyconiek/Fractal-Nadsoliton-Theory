import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3026_s1976_dissipation_time_order_reconciliation_no_new_live_frontier.py"
OUT = ROOT / "generated" / "p3026_s1976_dissipation_time_order_reconciliation_no_new_live_frontier.json"
MD = ROOT / "generated" / "p3026_s1976_dissipation_time_order_reconciliation_no_new_live_frontier.md"

class P3026DissipationTimeOrderReconciliationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P3026_DISSIPATION_TIME_ORDER_RECONCILIATION_NO_NEW_LIVE_FRONTIER_NO_CLOSURE")
        self.assertIsNotNone(self.payload["input_hashes"]["P3023"])
        self.assertIsNotNone(self.payload["input_hashes"]["P3024"])
        self.assertIsNotNone(self.payload["input_hashes"]["P3025"])

    def test_certificate_lattice(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["ledger_atom_count"], 3)
        self.assertEqual(cert["strict_source_closed_atoms"], 0)
        self.assertEqual(cert["closure_profile_count"], 8)
        self.assertEqual(cert["accepting_profile_count"], 1)
        self.assertFalse(cert["current_profile_accepts"])
        self.assertEqual(cert["new_live_frontier_count"], 0)

    def test_ledger_and_exports(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["object"], "DissipationTimeOrderLaneReconciliation_NoNewLiveFrontierCertificate")
        atoms = [row["atom"] for row in obj["ledger"]]
        self.assertEqual(atoms, ["directed_order_scaffold", "strict_chart_selector_source", "physical_tick_action_hamiltonian_unit"])
        self.assertTrue(all(not row["strict_source_closed"] for row in obj["ledger"]))
        self.assertTrue(all(not row["currently_supplied"] for row in obj["live_frontier_intake"]))
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3026/S1976", MD.read_text(encoding="utf-8"))
        self.assertIn("P3026/S1976", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3026/S1976", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3026", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
