import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2994_s1944_fourier_character_source_coupling_obstruction.py"
OUT = ROOT / "generated" / "p2994_s1944_fourier_character_source_coupling_obstruction.json"
MD = ROOT / "generated" / "p2994_s1944_fourier_character_source_coupling_obstruction.md"

class P2994FourierCharacterSourceCouplingObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2994_FOURIER_CHARACTER_SOURCE_COUPLING_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P2993"])

    def test_coupling_certificate(self):
        cert = self.payload["coupling_certificate"]
        self.assertEqual(cert["character_count"], 12)
        self.assertEqual(cert["named_source_atom_count"], 4)
        self.assertEqual(cert["coupling_test_count"], 48)
        self.assertTrue(cert["all_rows_have_exact_receivers"])
        self.assertEqual(cert["accepted_source_couplings"], [])
        self.assertEqual(cert["acceptance_matrix_rows"], 128)
        self.assertEqual(cert["accepted_rows"], 1)

    def test_obligations_and_rows(self):
        obj = self.payload["constructed_theoretical_objects"]
        obligations = {r["obligation"]: r["satisfied"] for r in obj["proof_obligation_rows"]}
        self.assertTrue(obligations["finite_character_atom_matrix_present"])
        self.assertTrue(obligations["exact_algebraic_receivers_present"])
        self.assertFalse(obligations["frequency_localizer_available"])
        self.assertFalse(obligations["strict_character_provenance_available"])
        self.assertFalse(obligations["atom_specific_coupling_theorem"])
        self.assertFalse(obligations["unit_bearing_coupling_coefficient"])
        self.assertFalse(obligations["nonproxy_export_available"])
        self.assertFalse(obligations["accepted_current_source_coupling"])
        rows = obj["coupling_witness"]["coupling_rows"]
        self.assertEqual(len(rows), 48)
        self.assertTrue(all(r["homomorphism_receiver_available"] for r in rows))
        self.assertFalse(any(r["accepted_source_coupling"] for r in rows))

    def test_docs_and_nonpromotion(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertIn("P2994/S1944", MD.read_text(encoding="utf-8"))
        self.assertIn("P2994/S1944", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2994/S1944", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2994", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
