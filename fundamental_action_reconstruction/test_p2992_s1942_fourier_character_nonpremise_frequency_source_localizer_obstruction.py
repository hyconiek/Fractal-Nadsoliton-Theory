import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2992_s1942_fourier_character_nonpremise_frequency_source_localizer_obstruction.py"
OUT = ROOT / "generated" / "p2992_s1942_fourier_character_nonpremise_frequency_source_localizer_obstruction.json"
MD = ROOT / "generated" / "p2992_s1942_fourier_character_nonpremise_frequency_source_localizer_obstruction.md"

class P2992FourierCharacterNonpremiseFrequencySourceLocalizerObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_input_and_theory_capsule(self):
        self.assertEqual(self.payload["status"], "P2992_FOURIER_CHARACTER_NONPREMISE_FREQUENCY_SOURCE_LOCALIZER_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P2991"])
        capsule = self.payload["theory_capsule"]
        self.assertIn("single primordial fractal information", capsule["ontology"])
        self.assertIn("sourced update/stationarity theorem", capsule["self_learning_reading"])

    def test_localizer_certificate(self):
        cert = self.payload["localizer_certificate"]
        self.assertEqual(cert["row_count"], 12)
        self.assertEqual(cert["unique_signature_rows"], [0, 6])
        self.assertEqual(cert["ambiguous_signature_classes"], [[1, 11], [2, 10], [3, 9], [4, 8], [5, 7]])
        self.assertEqual(cert["exported_frequency_source_localizers"], [])
        self.assertEqual(cert["acceptance_matrix_rows"], 128)
        self.assertEqual(cert["accepted_rows"], 1)

    def test_obligations_and_rows(self):
        obj = self.payload["constructed_theoretical_objects"]
        obligations = {r["obligation"]: r["satisfied"] for r in obj["proof_obligation_rows"]}
        self.assertTrue(obligations["finite_fourier_rows_present"])
        self.assertTrue(obligations["algebraic_frequency_signatures_constructed"])
        self.assertTrue(obligations["full_signature_has_unique_rows"])
        self.assertFalse(obligations["nonpremise_physical_frequency_sector"])
        self.assertFalse(obligations["strict_nadsoliton_character_source_map"])
        self.assertFalse(obligations["frequency_source_localizer_exported"])
        self.assertFalse(obligations["accepted_current_frequency_source_localizer"])
        rows = obj["localizer_witness"]["localizer_rows"]
        self.assertEqual(len(rows), 12)
        self.assertFalse(any(r["frequency_source_localizer_exported"] for r in rows))

    def test_docs_and_nonpromotion(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertIn("P2992/S1942", MD.read_text(encoding="utf-8"))
        self.assertIn("P2992/S1942", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2992/S1942", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2992", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
