import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3074_s2024_lyapunov_entropy_monotonicity_certificate.py"
OUT = ROOT / "generated" / "p3074_s2024_lyapunov_entropy_monotonicity_certificate.json"
MD = ROOT / "generated" / "p3074_s2024_lyapunov_entropy_monotonicity_certificate.md"

class P3074LyapunovEntropyMonotonicityCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input_hash(self):
        self.assertEqual(self.payload["status"], "P3074_INTERNAL_LYAPUNOV_MONOTONICITY_CERTIFICATE_NO_VARIATIONAL_PHYSICS_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P3073"])

    def test_matrix_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 3)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["p3073_accepted_intrinsic_bounded_scale_flow_rows"], 192)
        self.assertEqual(cert["accepted_p3073_flow_rows_reused"], 192)
        self.assertEqual(cert["step_sizes"], 3)
        self.assertEqual(cert["functionals_tested"], 3)
        self.assertEqual(cert["monotonicity_matrix_rows"], 1728)
        self.assertEqual(cert["accepted_internal_lyapunov_rows"], 1008)
        self.assertEqual(cert["variance_accepted_rows"], 528)
        self.assertEqual(cert["dirichlet_accepted_rows"], 480)
        self.assertEqual(cert["shell_control_monotone_rows"], 448)
        self.assertEqual(cert["shell_control_accepted_rows"], 0)
        self.assertEqual(cert["satisfied_proof_obligations"], 4)

    def test_functional_aggregates(self):
        rows = self.payload["constructed_theoretical_objects"]["functional_aggregate_certificate"]
        by_id = {row["functional"]: row for row in rows}
        self.assertEqual(by_id["variance_energy"]["accepted_internal_lyapunov_rows"], 528)
        self.assertEqual(by_id["quadratic_dirichlet_energy"]["accepted_internal_lyapunov_rows"], 480)
        self.assertEqual(by_id["shell_quadratic_energy"]["accepted_internal_lyapunov_rows"], 0)
        self.assertEqual(by_id["shell_quadratic_energy"]["monotone_nonincreasing_rows"], 448)

    def test_flags_recommendation_and_docs(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertTrue(self.payload["decision"]["positive_scoped_flags"]["internal_lyapunov_monotonicity_certificate_exported"])
        self.assertIn("bounded variational-source obstruction table", self.payload["decision"]["next_honest_step"])
        self.assertIn("P3074/S2024", MD.read_text(encoding="utf-8"))
        self.assertIn("P3074/S2024", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3074/S2024", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3074", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
