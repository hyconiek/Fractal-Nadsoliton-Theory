import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2980_s1930_z12_crt_idempotent_decomposition_source_candidate.py"
OUT = ROOT / "generated" / "p2980_s1930_z12_crt_idempotent_decomposition_source_candidate.json"
MD = ROOT / "generated" / "p2980_s1930_z12_crt_idempotent_decomposition_source_candidate.md"

class P2980Z12CRTIdempotentDecompositionSourceCandidateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2980_Z12_CRT_IDEMPOTENT_DECOMPOSITION_SOURCE_CANDIDATE_DEVELOPMENTAL_NO_STRICT_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P2979"])

    def test_idempotent_summary(self):
        summary = self.payload["idempotent_summary"]
        self.assertEqual(summary["idempotents"], [0, 1, 4, 9])
        self.assertEqual(summary["nontrivial_idempotents"], [4, 9])
        self.assertEqual(summary["orthogonal_completion_pair_count"], 1)
        pair = summary["orthogonal_completion_pairs"][0]
        self.assertEqual((pair["a"], pair["b"], pair["product_mod_12"], pair["sum_mod_12"]), (4, 9, 0, 1))
        self.assertEqual(summary["accepted_current_strict_source_objects"], [])
        self.assertEqual(summary["acceptance_matrix_rows"], 128)
        self.assertEqual(summary["accepted_rows"], 1)

    def test_obligations_and_rows(self):
        obligations = {r["obligation"]: r["satisfied"] for r in self.payload["constructed_theoretical_objects"]["proof_obligation_rows"]}
        self.assertTrue(obligations["new_outside_nilradical_lane"])
        self.assertTrue(obligations["finite_idempotent_witness"])
        self.assertTrue(obligations["orthogonal_projector_pair"])
        self.assertFalse(obligations["strict_nadsoliton_provenance_exported"])
        self.assertFalse(obligations["nonpremise_factor_semantics_exported"])
        self.assertFalse(obligations["couples_to_named_source_atom"])
        self.assertFalse(obligations["action_variational_installation"])
        self.assertFalse(obligations["accepted_current_strict_source_object"])
        rows = {r["candidate"]: r for r in self.payload["constructed_theoretical_objects"]["candidate_rows"]}
        self.assertTrue(rows["Z12_CRT_idempotent_decomposition_object"]["orthogonal_projector_pair"])
        self.assertFalse(rows["completed_strict_CRT_idempotent_source_schema"]["accepted_current_strict_source_object"])

    def test_docs_and_nonpromotion(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertIn("P2980/S1930", MD.read_text(encoding="utf-8"))
        self.assertIn("P2980/S1930", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2980/S1930", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2980", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
