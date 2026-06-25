import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3091_s2041_spectral_action_effective_action_obstruction_audit.py"
OUT = ROOT / "generated" / "p3091_s2041_spectral_action_effective_action_obstruction_audit.json"
MD = ROOT / "generated" / "p3091_s2041_spectral_action_effective_action_obstruction_audit.md"

class P3091SpectralActionEffectiveActionObstructionAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input_hash(self):
        self.assertEqual(self.payload["status"], "P3091_SPECTRAL_ACTION_EFFECTIVE_ACTION_GENERATING_FUNCTIONAL_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P3090"])

    def test_finite_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 3)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["p3090_accepted_nonimported_spectral_correlation_green_response_sources"], 0)
        self.assertEqual(cert["logdet_spectral_action_rows"], 3)
        self.assertEqual(cert["logdet_positive_eigenvalue_failures"], 0)
        self.assertEqual(cert["logdet_rows_with_action_unit"], 0)
        self.assertEqual(cert["source_coupled_quadratic_generator_rows"], 36)
        self.assertEqual(cert["quadratic_rows_with_unit_normalized_coupling"], 0)
        self.assertEqual(cert["quadratic_rows_with_empirical_response_generator"], 0)
        self.assertEqual(cert["coupling_scale_orbit_rows"], 4)
        self.assertEqual(cert["coupling_rows_with_unit_source"], 0)
        self.assertEqual(cert["finite_formal_variation_rows"], 12)
        self.assertEqual(cert["variation_rows_with_eom_source"], 0)
        self.assertEqual(cert["action_candidates"], 5)
        self.assertEqual(cert["candidate_gate_rows"], 30)
        self.assertEqual(cert["accepted_nonimported_spectral_action_effective_action_sources"], 0)
        self.assertEqual(cert["satisfied_proof_obligations"], 4)

    def test_witnesses_and_candidates(self):
        objs = self.payload["constructed_theoretical_objects"]
        self.assertTrue(all(row["positive_eigenvalues"] for row in objs["mass_regularized_logdet_spectral_action_rows"]))
        self.assertTrue(all(row["finite_source_coupled_quadratic_witness"] for row in objs["source_coupled_quadratic_generator_rows"]))
        self.assertTrue(all(row["quadratic_coupling_scaling_witness"] for row in objs["coupling_scale_orbit_rows"]))
        self.assertTrue(all(row["finite_formal_variation_witness"] for row in objs["finite_formal_variation_rows"]))
        aggs = objs["candidate_aggregate_certificate"]
        self.assertTrue(all(not row["accepted_nonimported_spectral_action_effective_action_source"] for row in aggs))
        by_id = {row["candidate"]: row for row in aggs}
        self.assertEqual(by_id["finite_z12_logdet_spectral_action"]["passed_gates"], 1)
        self.assertEqual(by_id["source_coupled_quadratic_form_JGJ"]["passed_gates"], 2)
        self.assertEqual(by_id["imported_empirical_effective_field_theory_template"]["passed_gates"], 5)

    def test_flags_recommendation_and_docs(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertTrue(self.payload["decision"]["positive_scoped_flags"]["source_coupled_quadratic_generator_computed"])
        self.assertIn("renormalization/scale-flow", self.payload["decision"]["next_honest_step"])
        self.assertIn("P3091/S2041", MD.read_text(encoding="utf-8"))
        self.assertIn("P3091/S2041", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3091/S2041", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3091", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
