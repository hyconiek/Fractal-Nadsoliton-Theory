import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3075_s2025_variational_source_obstruction_table.py"
OUT = ROOT / "generated" / "p3075_s2025_variational_source_obstruction_table.json"
MD = ROOT / "generated" / "p3075_s2025_variational_source_obstruction_table.md"

class P3075VariationalSourceObstructionTableTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input_hash(self):
        self.assertEqual(self.payload["status"], "P3075_LOCAL_DIRICHLET_VARIATIONAL_SOURCE_PARTIAL_EXPORT_MEAN_CENTERING_NONLOCAL_OBSTRUCTION")
        self.assertIsNotNone(self.payload["input_hashes"]["P3074"])

    def test_matrix_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 4)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["p3074_accepted_internal_lyapunov_rows"], 1008)
        self.assertEqual(cert["accepted_p3073_flow_rows_reused"], 192)
        self.assertEqual(cert["candidate_generators"], 5)
        self.assertEqual(cert["variational_source_matrix_rows"], 960)
        self.assertEqual(cert["exact_negative_gradient_match_rows"], 192)
        self.assertEqual(cert["accepted_local_variational_source_rows"], 96)
        self.assertEqual(cert["local_dirichlet_accepted_rows"], 96)
        self.assertEqual(cert["global_variance_exact_but_nonlocal_rows"], 96)
        self.assertEqual(cert["mean_centering_local_source_rows"], 0)
        self.assertEqual(cert["satisfied_proof_obligations"], 4)

    def test_generator_aggregates(self):
        rows = self.payload["constructed_theoretical_objects"]["generator_aggregate_certificate"]
        by_id = {row["candidate_generator"]: row for row in rows}
        self.assertEqual(by_id["local_dirichlet_generator"]["exact_match_rows"], 96)
        self.assertEqual(by_id["local_dirichlet_generator"]["accepted_local_variational_source_rows"], 96)
        self.assertEqual(by_id["global_variance_generator"]["exact_match_rows"], 96)
        self.assertEqual(by_id["global_variance_generator"]["accepted_local_variational_source_rows"], 0)
        self.assertEqual(by_id["linear_total_generator"]["total_preserving_source_rows"], 0)

    def test_flags_recommendation_and_docs(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertTrue(self.payload["decision"]["positive_scoped_flags"]["internal_local_dirichlet_gradient_source_for_laplacian_exported"])
        self.assertIn("continuum-limit/spectral-dispersion audit", self.payload["decision"]["next_honest_step"])
        self.assertIn("P3075/S2025", MD.read_text(encoding="utf-8"))
        self.assertIn("P3075/S2025", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3075/S2025", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3075", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
