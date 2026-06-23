import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3030_s1980_matter_spectral_field_localizer_obstruction.py"
OUT = ROOT / "generated" / "p3030_s1980_matter_spectral_field_localizer_obstruction.json"
MD = ROOT / "generated" / "p3030_s1980_matter_spectral_field_localizer_obstruction.md"

class P3030MatterSpectralFieldLocalizerTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3030_MATTER_SPECTRAL_FIELD_LOCALIZER_OBSTRUCTION_NO_MATTER_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P3029"])

    def test_orbit_degeneracy_and_no_localizer(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["dihedral_rows"], 24)
        self.assertEqual(cert["signature_preserving_rows"], 24)
        self.assertEqual(cert["translation_rows"], 12)
        self.assertEqual(cert["reflection_translation_rows"], 12)
        self.assertEqual(cert["localizer_candidates"], 3)
        self.assertEqual(cert["accepted_localizers"], 0)
        self.assertFalse(cert["field_representation_localizer_exported"])

    def test_obligations_and_exports(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["object"], "MatterSpectralFieldLocalizer_PhaseRetrievalOrbitObstructionMatrix")
        obligations = {row["obligation"]: row["satisfied"] for row in obj["proof_obligations"]}
        self.assertTrue(obligations["attacks_single_P3029_missing_atom"])
        self.assertTrue(obligations["uses_existing_P3029_generator"])
        self.assertFalse(obligations["translation_orbit_degeneracy_absent"])
        self.assertFalse(obligations["reflection_orbit_degeneracy_absent"])
        self.assertFalse(obligations["nonpremise_phase_recovery"])
        self.assertFalse(obligations["absolute_field_site_or_sector_localized"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3030/S1980", MD.read_text(encoding="utf-8"))
        self.assertIn("P3030/S1980", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3030/S1980", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3030", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
