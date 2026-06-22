import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3024_s1974_dissipation_chart_selector_source_obstruction.py"
OUT = ROOT / "generated" / "p3024_s1974_dissipation_chart_selector_source_obstruction.json"
MD = ROOT / "generated" / "p3024_s1974_dissipation_chart_selector_source_obstruction.md"

class P3024DissipationChartSelectorSourceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3024_DISSIPATION_CHART_SELECTOR_SOURCE_OBSTRUCTION_NO_CLOSURE")
        self.assertIsNotNone(self.payload["input_hashes"]["P3023"])

    def test_orbit_obstruction(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["u12_unit_count"], 4)
        self.assertEqual(cert["orbit_size"], 4)
        self.assertEqual(cert["stabilizer_size"], 1)
        self.assertEqual(cert["fixed_chart_rows"], 1)
        self.assertFalse(cert["accepted_as_strict_chart_selector_source"])

    def test_obligations_and_negative_exports(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["object"], "DissipationChainChartSelectorSource_AutOrbitObstructionMatrix")
        obligations = {row["obligation"]: row["satisfied"] for row in obj["proof_obligations"]}
        self.assertTrue(obligations["attacks_single_P3023_missing_theorem"])
        self.assertTrue(obligations["finite_U12_orbit_constructed"])
        self.assertFalse(obligations["nontrivial_stabilizer_for_chart"])
        self.assertFalse(obligations["U12_invariant_representative_exists"])
        self.assertFalse(obligations["endpoint_or_steepest_anchor_sources_chart"])
        self.assertFalse(obligations["strict_nonpremise_chart_selector_exported"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3024/S1974", MD.read_text(encoding="utf-8"))
        self.assertIn("P3024/S1974", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3024/S1974", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3024", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
