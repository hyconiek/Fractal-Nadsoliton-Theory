import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3090_s2040_spectral_correlation_green_function_response_obstruction_audit.py"
OUT = ROOT / "generated" / "p3090_s2040_spectral_correlation_green_function_response_obstruction_audit.json"
MD = ROOT / "generated" / "p3090_s2040_spectral_correlation_green_function_response_obstruction_audit.md"

class P3090SpectralCorrelationGreenFunctionResponseObstructionAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input_hash(self):
        self.assertEqual(self.payload["status"], "P3090_SPECTRAL_CORRELATION_GREEN_FUNCTION_RESPONSE_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P3089"])

    def test_finite_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 3)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["p3089_accepted_nonimported_born_rule_probability_readout_sources"], 0)
        self.assertEqual(cert["laplacian_spectrum_rows"], 12)
        self.assertEqual(cert["nonnegative_spectrum_failures"], 0)
        self.assertEqual(cert["zero_mode_rows"], 1)
        self.assertEqual(cert["pseudoinverse_green_kernel_rows"], 12)
        self.assertEqual(cert["green_rows_with_retarded_prescription"], 0)
        self.assertEqual(cert["mass_regularized_resolvent_rows"], 48)
        self.assertEqual(cert["resolvent_rows_with_unit_calibrated_spectral_density"], 0)
        self.assertEqual(cert["resolvent_rows_with_empirical_scattering_readout"], 0)
        self.assertEqual(cert["formal_iepsilon_regulator_rows"], 36)
        self.assertEqual(cert["iepsilon_rows_with_causal_retarded_source"], 0)
        self.assertEqual(cert["p3089_weighted_modal_correlation_rows"], 12)
        self.assertEqual(cert["modal_correlation_rows_with_response_law"], 0)
        self.assertEqual(cert["response_candidates"], 5)
        self.assertEqual(cert["candidate_gate_rows"], 30)
        self.assertEqual(cert["accepted_nonimported_spectral_correlation_green_response_sources"], 0)
        self.assertEqual(cert["satisfied_proof_obligations"], 4)

    def test_witnesses_and_candidates(self):
        objs = self.payload["constructed_theoretical_objects"]
        self.assertTrue(all(row["nonnegative_spectrum_witness"] for row in objs["laplacian_spectrum_rows"]))
        self.assertTrue(all(row["translation_invariant_kernel_witness"] for row in objs["pseudoinverse_green_kernel_rows"]))
        self.assertTrue(all(row["finite_response_kernel_witness"] for row in objs["mass_regularized_resolvent_rows"]))
        self.assertTrue(all(row["formal_bounded_regulator_witness"] for row in objs["formal_iepsilon_regulator_rows"]))
        self.assertTrue(all(row["two_point_correlation_witness"] for row in objs["p3089_weighted_modal_correlation_rows"]))
        aggs = objs["candidate_aggregate_certificate"]
        self.assertTrue(all(not row["accepted_nonimported_spectral_correlation_green_response_source"] for row in aggs))
        by_id = {row["candidate"]: row for row in aggs}
        self.assertEqual(by_id["z12_laplacian_pseudoinverse_green_kernel"]["passed_gates"], 2)
        self.assertEqual(by_id["p3089_probability_weighted_modal_correlator"]["passed_gates"], 1)
        self.assertEqual(by_id["imported_scattering_spectral_density_template"]["passed_gates"], 5)

    def test_flags_recommendation_and_docs(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertTrue(self.payload["decision"]["positive_scoped_flags"]["pseudoinverse_green_kernel_computed"])
        self.assertIn("effective-action", self.payload["decision"]["next_honest_step"])
        self.assertIn("P3090/S2040", MD.read_text(encoding="utf-8"))
        self.assertIn("P3090/S2040", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3090/S2040", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3090", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
