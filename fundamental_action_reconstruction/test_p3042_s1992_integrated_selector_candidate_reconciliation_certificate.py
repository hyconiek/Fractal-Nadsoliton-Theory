import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3042_s1992_integrated_selector_candidate_reconciliation_certificate.py"
OUT = ROOT / "generated" / "p3042_s1992_integrated_selector_candidate_reconciliation_certificate.json"
MD = ROOT / "generated" / "p3042_s1992_integrated_selector_candidate_reconciliation_certificate.md"

class P3042IntegratedSelectorCandidateReconciliationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(cls_status := self.payload["status"], "P3042_INTEGRATED_SELECTOR_CANDIDATE_RECONCILIATION_NO_SELECTOR_EXPORT")
        self.assertTrue(cls_status.endswith("NO_SELECTOR_EXPORT"))
        for key in ["P3038", "P3039", "P3040", "P3041"]:
            self.assertIsNotNone(self.payload["input_hashes"][key])

    def test_lattice_certificate(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["atom_count"], 5)
        self.assertEqual(cert["closed_atoms"], 1)
        self.assertEqual(cert["open_required_source_atoms"], 4)
        self.assertEqual(cert["profile_count"], 32)
        self.assertEqual(cert["accepted_counterfactual_profiles"], 1)
        self.assertFalse(cert["current_profile_accepted"])
        self.assertFalse(cert["strict_selector_mechanism_exported"])

    def test_obligations_and_exports(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["object"], "IntegratedSelectorCandidate_ReconciliationNoExportCertificate")
        obligations = {row["obligation"]: row["satisfied"] for row in obj["proof_obligations"]}
        self.assertTrue(obligations["p3038_integrated_operator_accounted"])
        self.assertTrue(obligations["p3039_chiral_source_audited"])
        self.assertTrue(obligations["p3040_path_source_audited"])
        self.assertTrue(obligations["p3041_unit_readout_audited"])
        self.assertFalse(obligations["all_required_source_atoms_closed"])
        self.assertFalse(obligations["current_profile_accepts_selector_export"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3042/S1992", MD.read_text(encoding="utf-8"))
        self.assertIn("P3042/S1992", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3042/S1992", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3042", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
