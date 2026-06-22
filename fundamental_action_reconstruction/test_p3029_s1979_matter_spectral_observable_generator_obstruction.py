import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3029_s1979_matter_spectral_observable_generator_obstruction.py"
OUT = ROOT / "generated" / "p3029_s1979_matter_spectral_observable_generator_obstruction.json"
MD = ROOT / "generated" / "p3029_s1979_matter_spectral_observable_generator_obstruction.md"

class P3029MatterSpectralObservableGeneratorTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3029_MATTER_SPECTRAL_OBSERVABLE_GENERATOR_OBSTRUCTION_NO_CLASSICAL_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P3028"])

    def test_generator_positive_but_no_matter_export(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["signature_length"], 12)
        self.assertEqual(cert["unit_rows"], 4)
        self.assertEqual(cert["unit_invariant_rows"], 4)
        self.assertTrue(cert["observer_independent_generator_accepted"])
        self.assertFalse(cert["matter_sector_export_accepted"])

    def test_obligations_and_exports(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["object"], "MatterSpectralObservableGenerator_DFTMagnitudeInvariantObstructionMatrix")
        obligations = {row["obligation"]: row["satisfied"] for row in obj["proof_obligations"]}
        self.assertTrue(obligations["attacks_single_P3028_foundation_atom"])
        self.assertTrue(obligations["explicit_domain_codomain"])
        self.assertTrue(obligations["observer_independent_formula"])
        self.assertTrue(obligations["U12_unit_invariant_signature"])
        self.assertFalse(obligations["field_representation_localizer"])
        self.assertFalse(obligations["mass_coupling_provenance"])
        self.assertFalse(obligations["selector_or_sector_source"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3029/S1979", MD.read_text(encoding="utf-8"))
        self.assertIn("P3029/S1979", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3029/S1979", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3029", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
