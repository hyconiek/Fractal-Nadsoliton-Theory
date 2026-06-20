import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2988_s1938_annihilator_lattice_nonpremise_source_localizer_obstruction.py"
OUT = ROOT / "generated" / "p2988_s1938_annihilator_lattice_nonpremise_source_localizer_obstruction.json"
MD = ROOT / "generated" / "p2988_s1938_annihilator_lattice_nonpremise_source_localizer_obstruction.md"

class P2988AnnihilatorLatticeNonpremiseSourceLocalizerObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2988_ANNIHILATOR_LATTICE_NONPREMISE_SOURCE_LOCALIZER_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P2987"])

    def test_localizer_certificate(self):
        cert = self.payload["localizer_certificate"]
        self.assertEqual(cert["row_count"], 6)
        self.assertTrue(cert["all_rows_algebraically_distinguished"])
        self.assertEqual(cert["exported_nonpremise_source_localizers"], [])
        self.assertEqual(cert["acceptance_matrix_rows"], 64)
        self.assertEqual(cert["accepted_rows"], 1)

    def test_obligations_and_rows(self):
        obj = self.payload["constructed_theoretical_objects"]
        obligations = {r["obligation"]: r["satisfied"] for r in obj["proof_obligation_rows"]}
        self.assertTrue(obligations["finite_annihilator_rows_present"])
        self.assertTrue(obligations["algebraic_signature_distinguishes_rows"])
        self.assertFalse(obligations["nonpremise_physical_sector"])
        self.assertFalse(obligations["strict_nadsoliton_row_selector"])
        self.assertFalse(obligations["source_localizer_exported"])
        self.assertFalse(obligations["accepted_current_localizer"])
        rows = obj["localizer_witness"]["localizer_rows"]
        self.assertEqual(len(rows), 6)
        self.assertTrue(all(r["row_distinguished_by_signature"] for r in rows))
        self.assertFalse(any(r["source_localizer_exported"] for r in rows))

    def test_docs_and_nonpromotion(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertIn("P2988/S1938", MD.read_text(encoding="utf-8"))
        self.assertIn("P2988/S1938", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2988/S1938", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2988", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
