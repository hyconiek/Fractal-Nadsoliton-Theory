import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3031_s1981_matter_spectral_mass_coupling_provenance_obstruction.py"
OUT = ROOT / "generated" / "p3031_s1981_matter_spectral_mass_coupling_provenance_obstruction.json"
MD = ROOT / "generated" / "p3031_s1981_matter_spectral_mass_coupling_provenance_obstruction.md"

class P3031MatterSpectralMassCouplingProvenanceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3031_MATTER_SPECTRAL_MASS_COUPLING_PROVENANCE_OBSTRUCTION_NO_MATTER_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P3030"])

    def test_rescaling_and_candidates(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["rescaling_test"]["K_scale_factor"], 3.0)
        self.assertEqual(obj["rescaling_test"]["raw_signature_l1_scale_ratio"], 3.0)
        self.assertTrue(obj["rescaling_test"]["normalized_signature_preserved"])
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["candidate_rows"], 4)
        self.assertEqual(cert["accepted_candidates"], 0)
        self.assertGreaterEqual(cert["raw_scale_covariant_rows"], 2)
        self.assertEqual(cert["dimensionless_unitless_rows"], 2)
        self.assertFalse(cert["mass_coupling_provenance_exported"])

    def test_obligations_and_exports(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["object"], "MatterSpectralMassCouplingProvenance_RescalingNormalizationObstructionMatrix")
        obligations = {row["obligation"]: row["satisfied"] for row in obj["proof_obligations"]}
        self.assertTrue(obligations["attacks_single_P3029_missing_atom"])
        self.assertTrue(obligations["uses_existing_P3029_generator"])
        self.assertFalse(obligations["field_representation_localizer_available"])
        self.assertFalse(obligations["absolute_mass_or_coupling_unit_source"])
        self.assertFalse(obligations["target_independent_mass_coupling_map"])
        self.assertFalse(obligations["unit_bearing_action_eom_insertion"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3031/S1981", MD.read_text(encoding="utf-8"))
        self.assertIn("P3031/S1981", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3031/S1981", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3031", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
