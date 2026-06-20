import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2993_s1943_fourier_character_strict_provenance_obstruction.py"
OUT = ROOT / "generated" / "p2993_s1943_fourier_character_strict_provenance_obstruction.json"
MD = ROOT / "generated" / "p2993_s1943_fourier_character_strict_provenance_obstruction.md"

class P2993FourierCharacterStrictProvenanceObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2993_FOURIER_CHARACTER_STRICT_PROVENANCE_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P2992"])

    def test_provenance_certificate(self):
        cert = self.payload["provenance_certificate"]
        self.assertEqual(cert["row_count"], 12)
        self.assertTrue(cert["all_rows_exact_additive_characters"])
        self.assertTrue(cert["p2991_orthogonality_retained"])
        self.assertEqual(cert["accepted_strict_character_provenance_rows"], [])
        self.assertEqual(cert["acceptance_matrix_rows"], 128)
        self.assertEqual(cert["accepted_rows"], 1)

    def test_obligations_and_rows(self):
        obj = self.payload["constructed_theoretical_objects"]
        obligations = {r["obligation"]: r["satisfied"] for r in obj["proof_obligation_rows"]}
        self.assertTrue(obligations["finite_character_rows_present"])
        self.assertTrue(obligations["exact_additive_homomorphism_rows"])
        self.assertTrue(obligations["orthogonality_certificate_retained"])
        self.assertFalse(obligations["frequency_localizer_available"])
        self.assertFalse(obligations["strict_nadsoliton_source_map_exported"])
        self.assertFalse(obligations["nonpremise_internal_character_provenance"])
        self.assertFalse(obligations["accepted_current_strict_character_provenance"])
        rows = obj["provenance_witness"]["provenance_rows"]
        self.assertEqual(len(rows), 12)
        self.assertTrue(all(r["homomorphism_defect_count"] == 0 for r in rows))
        self.assertFalse(any(r["accepted_strict_character_provenance"] for r in rows))

    def test_docs_and_nonpromotion(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertIn("P2993/S1943", MD.read_text(encoding="utf-8"))
        self.assertIn("P2993/S1943", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2993/S1943", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2993", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
