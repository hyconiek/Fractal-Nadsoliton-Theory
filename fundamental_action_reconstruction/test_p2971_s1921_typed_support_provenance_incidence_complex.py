import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2971_s1921_typed_support_provenance_incidence_complex.py"
OUT = ROOT / "generated" / "p2971_s1921_typed_support_provenance_incidence_complex.json"
MD = ROOT / "generated" / "p2971_s1921_typed_support_provenance_incidence_complex.md"

class P2971TypedSupportProvenanceIncidenceComplexTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2971_TYPED_SUPPORT_PROVENANCE_INCIDENCE_COMPLEX_DEVELOPMENTAL_NO_STRICT_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P2970"])

    def test_structural_certificate(self):
        cert = self.payload["structural_certificate"]
        self.assertEqual(cert["aggregate_vector"], [1, 2, 2, 2, 2])
        self.assertEqual(cert["aggregate_sum"], 9)
        self.assertEqual(cert["primitive_mean"], "9/5")
        self.assertEqual(cert["component_weight_automorphisms"], 6)
        self.assertEqual(cert["provenance_preserving_automorphisms"], 1)
        self.assertEqual(cert["acceptance_matrix_rows"], 128)
        self.assertEqual(cert["strict_rows"], 1)

    def test_obligations_and_rows(self):
        obligations = {r["obligation"]: r["satisfied"] for r in self.payload["constructed_theoretical_objects"]["proof_obligation_rows"]}
        self.assertTrue(obligations["new_typed_structural_object_constructed"])
        self.assertTrue(obligations["outside_ratio_package_arithmetic"])
        self.assertTrue(obligations["KC_mismatch_preserved"])
        self.assertTrue(obligations["provenance_noncollapse_checked"])
        self.assertFalse(obligations["strict_source_localizer_exported"])
        self.assertFalse(obligations["unit_bearing_coupling_exported"])
        self.assertFalse(obligations["nonproxy_variational_chain_rule_exported"])
        rows = {r["object"]: r for r in self.payload["constructed_theoretical_objects"]["structural_rows"]}
        self.assertTrue(rows["typed_support_provenance_incidence_complex"]["outside_ratio_package_arithmetic"])
        self.assertFalse(rows["typed_support_provenance_incidence_complex"]["accepted_current_strict_object"])

    def test_docs_and_nonpromotion(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertIn("P2971/S1921", MD.read_text(encoding="utf-8"))
        self.assertIn("P2971/S1921", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2971/S1921", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2971", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
