import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3035_s1985_z12_directional_branch_selector_source_obstruction.py"
OUT = ROOT / "generated" / "p3035_s1985_z12_directional_branch_selector_source_obstruction.json"
MD = ROOT / "generated" / "p3035_s1985_z12_directional_branch_selector_source_obstruction.md"

class P3035Z12DirectionalBranchSelectorSourceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P3035_Z12_DIRECTIONAL_BRANCH_SELECTOR_SOURCE_OBSTRUCTION_NO_SELECTOR_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P3028"])
        self.assertIsNotNone(self.payload["input_hashes"]["P3034"])

    def test_finite_certificate(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["candidate_rows"], 4)
        self.assertEqual(cert["orientation_separating_rows"], 0)
        self.assertEqual(cert["accepted_branch_source_rows"], 0)
        self.assertEqual(cert["inversion_units_exchange_branches"], 2)
        self.assertFalse(cert["strict_selector_branch_source_exported"])

    def test_obligations_and_exports(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["object"], "Z12DirectionalBranchSelectorSource_ObstructionMatrix")
        obligations = {row["obligation"]: row["satisfied"] for row in obj["proof_obligations"]}
        self.assertTrue(obligations["attacks_single_P3028_foundation_atom"])
        self.assertTrue(obligations["explicit_orientation_torsor"])
        self.assertTrue(obligations["finite_candidate_receivers_computable"])
        self.assertFalse(obligations["orientation_score_separation"])
        self.assertFalse(obligations["aut_z12_compatible_nonpremise_source"])
        self.assertFalse(obligations["chart_independent_branch_localizer"])
        self.assertFalse(obligations["coupling_to_classical_readout_rows"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3035/S1985", MD.read_text(encoding="utf-8"))
        self.assertIn("P3035/S1985", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3035/S1985", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3035", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
