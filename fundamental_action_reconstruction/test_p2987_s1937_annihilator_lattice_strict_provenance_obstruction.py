import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2987_s1937_annihilator_lattice_strict_provenance_obstruction.py"
OUT = ROOT / "generated" / "p2987_s1937_annihilator_lattice_strict_provenance_obstruction.json"
MD = ROOT / "generated" / "p2987_s1937_annihilator_lattice_strict_provenance_obstruction.md"

class P2987AnnihilatorLatticeStrictProvenanceObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2987_ANNIHILATOR_LATTICE_STRICT_PROVENANCE_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P2986"])

    def test_provenance_certificate(self):
        cert = self.payload["provenance_certificate"]
        self.assertEqual(cert["ideal_count"], 6)
        self.assertTrue(cert["finite_annihilator_rows_exact"])
        self.assertTrue(cert["order_reversing_checks_pass"])
        self.assertFalse(cert["proper_rows_translation_preserved"])
        self.assertEqual(cert["accepted_current_provenance_rows"], [])
        self.assertEqual(cert["acceptance_matrix_rows"], 64)
        self.assertEqual(cert["accepted_rows"], 1)

    def test_obligations_and_rows(self):
        obj = self.payload["constructed_theoretical_objects"]
        obligations = {r["obligation"]: r["satisfied"] for r in obj["proof_obligation_rows"]}
        self.assertTrue(obligations["finite_annihilator_lattice_recomputed"])
        self.assertTrue(obligations["order_reversing_lattice_retained"])
        self.assertFalse(obligations["proper_ideal_translation_invariance"])
        self.assertFalse(obligations["strict_nadsoliton_source_map"])
        self.assertFalse(obligations["nonpremise_row_localizer"])
        self.assertFalse(obligations["accepted_strict_provenance"])
        rows = obj["strict_provenance_witness"]["provenance_rows"]
        self.assertEqual(len(rows), 6)
        self.assertTrue(all(r["finite_annihilator_exact"] for r in rows))
        self.assertFalse(any(r["accepted_strict_provenance"] for r in rows))

    def test_docs_and_nonpromotion(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertIn("P2987/S1937", MD.read_text(encoding="utf-8"))
        self.assertIn("P2987/S1937", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2987/S1937", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2987", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
