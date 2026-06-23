import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3036_s1986_external_physical_unit_source_scale_orbit_obstruction.py"
OUT = ROOT / "generated" / "p3036_s1986_external_physical_unit_source_scale_orbit_obstruction.json"
MD = ROOT / "generated" / "p3036_s1986_external_physical_unit_source_scale_orbit_obstruction.md"

class P3036ExternalPhysicalUnitSourceScaleOrbitTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P3036_EXTERNAL_PHYSICAL_UNIT_SOURCE_SCALE_ORBIT_OBSTRUCTION_NO_UNIT_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P3028"])
        self.assertIsNotNone(self.payload["input_hashes"]["P3035"])

    def test_finite_certificate(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["candidate_rows"], 4)
        self.assertEqual(cert["dimensionless_representative_rows"], 4)
        self.assertEqual(cert["scale_orbit_invariant_rows"], 0)
        self.assertEqual(cert["accepted_external_unit_source_rows"], 0)
        self.assertFalse(cert["external_physical_unit_source_exported"])

    def test_obligations_and_exports(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["object"], "ExternalPhysicalUnitSource_ScaleOrbitObstructionMatrix")
        obligations = {row["obligation"]: row["satisfied"] for row in obj["proof_obligations"]}
        self.assertTrue(obligations["attacks_single_P3028_foundation_atom"])
        self.assertTrue(obligations["explicit_scale_torsor"])
        self.assertTrue(obligations["finite_candidate_unit_receivers_computable"])
        self.assertFalse(obligations["dimensionful_external_unit_export"])
        self.assertFalse(obligations["amplitude_scale_orbit_quotient"])
        self.assertFalse(obligations["label_length_scale_orbit_quotient"])
        self.assertFalse(obligations["reference_cell_and_bit_to_unit_map"])
        self.assertFalse(obligations["coupling_to_unit_bearing_action_or_readout"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3036/S1986", MD.read_text(encoding="utf-8"))
        self.assertIn("P3036/S1986", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3036/S1986", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3036", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
