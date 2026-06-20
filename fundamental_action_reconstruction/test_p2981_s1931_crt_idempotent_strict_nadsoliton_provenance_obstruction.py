import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2981_s1931_crt_idempotent_strict_nadsoliton_provenance_obstruction.py"
OUT = ROOT / "generated" / "p2981_s1931_crt_idempotent_strict_nadsoliton_provenance_obstruction.json"
MD = ROOT / "generated" / "p2981_s1931_crt_idempotent_strict_nadsoliton_provenance_obstruction.md"

class P2981CRTIdempotentStrictNadsolitonProvenanceObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2981_CRT_IDEMPOTENT_STRICT_NADSOLITON_PROVENANCE_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P2980"])

    def test_provenance_certificate(self):
        cert = self.payload["provenance_certificate"]
        self.assertEqual(cert["idempotents"], [0, 1, 4, 9])
        self.assertTrue(cert["unique_nontrivial_completion_pair"])
        self.assertTrue(cert["factor_cardinality_distinguishes_labels"])
        self.assertEqual(cert["orthogonal_completion_pairs"][0]["a"], 4)
        self.assertEqual(cert["orthogonal_completion_pairs"][0]["b"], 9)
        self.assertEqual(cert["accepted_current_strict_provenance_theorems"], [])
        self.assertEqual(cert["acceptance_matrix_rows"], 64)
        self.assertEqual(cert["accepted_rows"], 1)

    def test_obligations_and_rows(self):
        obligations = {r["obligation"]: r["satisfied"] for r in self.payload["constructed_theoretical_objects"]["proof_obligation_rows"]}
        self.assertTrue(obligations["finite_projector_uniqueness"])
        self.assertTrue(obligations["factor_label_distinction"])
        self.assertFalse(obligations["strict_nadsoliton_source_map_exported"])
        self.assertFalse(obligations["nonpremise_factor_semantics_exported"])
        self.assertFalse(obligations["couples_to_named_source_atom"])
        self.assertFalse(obligations["accepted_strict_provenance_theorem"])
        rows = {r["candidate"]: r for r in self.payload["constructed_theoretical_objects"]["candidate_rows"]}
        self.assertTrue(rows["finite_CRT_projector_uniqueness"]["finite_projector_uniqueness"])
        self.assertFalse(rows["completed_strict_CRT_projector_provenance_schema"]["accepted_strict_provenance_theorem"])

    def test_docs_and_nonpromotion(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertIn("P2981/S1931", MD.read_text(encoding="utf-8"))
        self.assertIn("P2981/S1931", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2981/S1931", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2981", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
