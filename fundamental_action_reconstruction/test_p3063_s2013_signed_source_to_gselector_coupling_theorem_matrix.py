import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3063_s2013_signed_source_to_gselector_coupling_theorem_matrix.py"
OUT = ROOT / "generated" / "p3063_s2013_signed_source_to_gselector_coupling_theorem_matrix.json"
MD = ROOT / "generated" / "p3063_s2013_signed_source_to_gselector_coupling_theorem_matrix.md"

class P3063SignedSourceToGSelectorCouplingTheoremMatrixTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3063_SIGNED_SOURCE_TO_GSELECTOR_COUPLING_THEOREM_MATRIX_NO_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P3062"])

    def test_lay_analogy_boundary(self):
        boundary = self.payload["constructed_theoretical_objects"]["lay_analogy_boundary"]
        self.assertIn("operational", boundary["kernel_as_laws"])
        self.assertIn("weak", boundary["selector_as_start"])
        self.assertIn("nadsoliton", boundary["ontology_correction"])
        flags = self.payload["decision"]["negative_export_flags"]
        self.assertFalse(flags["simulation_start_exported"])
        self.assertFalse(flags["neural_network_universe_start_exported"])

    def test_finite_coupling_matrix(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 3)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["signed_source_rows"], 4)
        self.assertEqual(cert["coupling_polarities_per_source"], 2)
        self.assertEqual(cert["coupling_theorem_rows"], 8)
        self.assertEqual(cert["rows_with_coupling_map_exists"], 8)
        self.assertEqual(cert["rows_with_unique_polarity_selected_by_theory"], 0)
        self.assertEqual(cert["accepted_coupling_theorems"], 0)
        self.assertEqual(cert["satisfied_proof_obligations"], 3)

    def test_theorem_rows_and_negative_exports(self):
        obj = self.payload["constructed_theoretical_objects"]
        normal = obj["coupling_theorem_normal_form"]
        self.assertEqual(normal["object"], "ExplicitSignedSourceToGSelectorCouplingTheoremMatrix")
        self.assertEqual(len(normal["criteria"]), 6)
        rows = obj["theorem_rows"]
        self.assertEqual(len(rows), 8)
        self.assertTrue(all(row["criteria"]["coupling_map_exists"] for row in rows))
        self.assertTrue(all(not row["accepted_coupling_theorem"] for row in rows))
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3063/S2013", MD.read_text(encoding="utf-8"))
        self.assertIn("P3063/S2013", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3063/S2013", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3063", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
