import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2977_s1927_nilradical_strict_nadsoliton_provenance_obstruction.py"
OUT = ROOT / "generated" / "p2977_s1927_nilradical_strict_nadsoliton_provenance_obstruction.json"
MD = ROOT / "generated" / "p2977_s1927_nilradical_strict_nadsoliton_provenance_obstruction.md"

class P2977NilradicalStrictNadsolitonProvenanceObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2977_NILRADICAL_STRICT_NADSOLITON_PROVENANCE_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P2976"])

    def test_provenance_certificate(self):
        cert = self.payload["provenance_certificate"]
        self.assertEqual(cert["candidate_count"], 4)
        self.assertEqual(cert["translated_antipodal_coset_count"], 6)
        self.assertEqual(cert["translation_stabilizer"], [0, 6])
        self.assertTrue(cert["all_units_fix_nilradical_set"])
        self.assertEqual(cert["acceptance_matrix_rows"], 64)
        self.assertEqual(cert["accepted_rows"], 1)
        self.assertEqual(cert["accepted_current_strict_provenance_theorems"], [])

    def test_obligations_and_candidate_rows(self):
        obligations = {r["obligation"]: r["satisfied"] for r in self.payload["constructed_theoretical_objects"]["proof_obligation_rows"]}
        self.assertTrue(obligations["finite_ring_canonicity"])
        self.assertTrue(obligations["unit_invariance"])
        self.assertFalse(obligations["translation_orbit_degeneracy_removed"])
        self.assertFalse(obligations["strict_nadsoliton_source_map_exported"])
        self.assertFalse(obligations["nonpremise_zero_section_source_exported"])
        self.assertFalse(obligations["couples_to_named_missing_source_atom"])
        self.assertFalse(obligations["accepted_strict_provenance_theorem"])
        rows = {r["candidate"]: r for r in self.payload["constructed_theoretical_objects"]["candidate_rows"]}
        self.assertTrue(rows["ring_canonical_nilradical_provenance"]["finite_ring_canonicity"])
        self.assertFalse(rows["completed_strict_nilradical_provenance_schema"]["accepted_strict_provenance_theorem"])

    def test_docs_and_nonpromotion(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertIn("P2977/S1927", MD.read_text(encoding="utf-8"))
        self.assertIn("P2977/S1927", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2977/S1927", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2977", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
