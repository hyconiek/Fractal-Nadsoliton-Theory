import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3032_s1982_matter_spectral_lane_reconciliation_no_new_live_frontier.py"
OUT = ROOT / "generated" / "p3032_s1982_matter_spectral_lane_reconciliation_no_new_live_frontier.json"
MD = ROOT / "generated" / "p3032_s1982_matter_spectral_lane_reconciliation_no_new_live_frontier.md"

class P3032MatterSpectralLaneReconciliationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P3032_MATTER_SPECTRAL_LANE_RECONCILIATION_NO_NEW_LIVE_FRONTIER")
        self.assertIsNotNone(self.payload["input_hashes"]["P3029"])
        self.assertIsNotNone(self.payload["input_hashes"]["P3030"])
        self.assertIsNotNone(self.payload["input_hashes"]["P3031"])

    def test_ledger_lattice_and_frontier(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["ledger_atoms"], 4)
        self.assertEqual(cert["closed_atoms"], 1)
        self.assertEqual(cert["closure_profiles"], 16)
        self.assertEqual(cert["accepting_profiles"], 1)
        self.assertFalse(cert["current_profile_accepts_matter_sector"])
        self.assertEqual(cert["new_live_frontier_count"], 0)

    def test_atom_status_and_negative_exports(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["object"], "MatterSpectralLaneReconciliation_NoNewLiveFrontierCertificate")
        atom_status = {row["atom"]: row["strict_source_closed"] for row in obj["atom_ledger"]}
        self.assertTrue(atom_status["observer_independent_observable_generator"])
        self.assertFalse(atom_status["field_representation_localizer"])
        self.assertFalse(atom_status["mass_coupling_provenance"])
        self.assertFalse(atom_status["selector_sector_unit_action_insertion"])
        for row in obj["live_frontier_intake"]:
            self.assertFalse(row["currently_supplied"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3032/S1982", MD.read_text(encoding="utf-8"))
        self.assertIn("P3032/S1982", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3032/S1982", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3032", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
