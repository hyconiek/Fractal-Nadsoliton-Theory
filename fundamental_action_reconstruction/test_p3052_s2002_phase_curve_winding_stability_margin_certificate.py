import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3052_s2002_phase_curve_winding_stability_margin_certificate.py"
OUT = ROOT / "generated" / "p3052_s2002_phase_curve_winding_stability_margin_certificate.json"
MD = ROOT / "generated" / "p3052_s2002_phase_curve_winding_stability_margin_certificate.md"

class P3052PhaseCurveWindingStabilityMarginCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3052_PHASE_CURVE_WINDING_STABILITY_MARGIN_CERTIFICATE_BOUNDED_NO_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P3051"])

    def test_finite_certificate(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["base_integer_winding"], 1)
        self.assertEqual(cert["clearance_rows"], 12)
        self.assertEqual(cert["positive_clearance_rows"], 12)
        self.assertGreater(cert["minimum_centroid_edge_clearance"], 0)
        self.assertEqual(cert["perturbation_rows"], 40)
        self.assertEqual(cert["winding_preserved_perturbation_rows"], 40)
        self.assertEqual(cert["aut_signed_winding_rows"], 4)
        self.assertEqual(cert["aut_signed_winding_verified_rows"], 4)
        self.assertEqual(cert["source_acceptance_criteria"], 8)
        self.assertEqual(cert["satisfied_source_acceptance_criteria"], 4)
        self.assertFalse(cert["strict_winding_source_theorem_exported"])

    def test_obligations_and_exports(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["object"], "PhaseCurveWinding_StabilityMarginCertificate")
        obligations = {row["obligation"]: row["satisfied"] for row in obj["proof_obligations"]}
        self.assertTrue(obligations["p3051_winding_stability_targeted"])
        self.assertTrue(obligations["positive_clearance_certificate"])
        self.assertTrue(obligations["finite_perturbation_witness_table"])
        self.assertTrue(obligations["aut_signed_winding_boundary"])
        self.assertFalse(obligations["strict_orientation_source_theorem"])
        self.assertFalse(obligations["selector_ltotal_bridge_toe_installation"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3052/S2002", MD.read_text(encoding="utf-8"))
        self.assertIn("P3052/S2002", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3052/S2002", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3052", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
