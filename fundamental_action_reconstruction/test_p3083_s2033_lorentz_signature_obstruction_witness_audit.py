import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3083_s2033_lorentz_signature_obstruction_witness_audit.py"
OUT = ROOT / "generated" / "p3083_s2033_lorentz_signature_obstruction_witness_audit.json"
MD = ROOT / "generated" / "p3083_s2033_lorentz_signature_obstruction_witness_audit.md"

class P3083LorentzSignatureObstructionWitnessAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input_hash(self):
        self.assertEqual(self.payload["status"], "P3083_LORENTZ_SIGNATURE_OBSTRUCTION_WITNESS_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P3082"])

    def test_finite_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 3)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["p3082_accepted_nonimported_continuum_limit_functors"], 0)
        self.assertEqual(cert["quadratic_form_signature_rows"], 5)
        self.assertEqual(cert["internal_signature_rows"], 3)
        self.assertEqual(cert["internal_indefinite_signature_rows"], 0)
        self.assertEqual(cert["binary_profile_rows"], 4096)
        self.assertEqual(cert["binary_rows_with_nonnegative_dirichlet_energy"], 4096)
        self.assertEqual(cert["binary_rows_separating_time_axis"], 0)
        self.assertEqual(cert["signature_candidates"], 5)
        self.assertEqual(cert["candidate_gate_rows"], 30)
        self.assertEqual(cert["accepted_nonimported_lorentz_signature_sources"], 0)
        self.assertEqual(cert["satisfied_proof_obligations"], 4)

    def test_signature_rows_and_candidate_aggregates(self):
        rows = {row["form_id"]: row for row in self.payload["constructed_theoretical_objects"]["quadratic_form_signature_rows"]}
        self.assertEqual(rows["z12_laplacian_L"]["positive_count"], 11)
        self.assertEqual(rows["z12_laplacian_L"]["zero_count"], 1)
        self.assertFalse(rows["z12_laplacian_L"]["indefinite"])
        self.assertTrue(rows["formal_hyperbolic_block_minus_dt2_plus_L"]["indefinite"])
        aggs = self.payload["constructed_theoretical_objects"]["candidate_aggregate_certificate"]
        self.assertTrue(all(not row["accepted_nonimported_lorentz_signature_source"] for row in aggs))
        by_id = {row["candidate"]: row for row in aggs}
        self.assertEqual(by_id["z12_dirichlet_euclidean_form"]["passed_gates"], 1)
        self.assertEqual(by_id["imported_minkowski_metric_template"]["passed_gates"], 5)

    def test_flags_recommendation_and_docs(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertTrue(self.payload["decision"]["positive_scoped_flags"]["internal_quadratic_signatures_computed"])
        self.assertIn("gauge-representation", self.payload["decision"]["next_honest_step"])
        self.assertIn("P3083/S2033", MD.read_text(encoding="utf-8"))
        self.assertIn("P3083/S2033", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3083/S2033", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3083", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
