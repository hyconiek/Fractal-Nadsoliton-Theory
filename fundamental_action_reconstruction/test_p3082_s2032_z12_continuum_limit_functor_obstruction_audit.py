import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3082_s2032_z12_continuum_limit_functor_obstruction_audit.py"
OUT = ROOT / "generated" / "p3082_s2032_z12_continuum_limit_functor_obstruction_audit.json"
MD = ROOT / "generated" / "p3082_s2032_z12_continuum_limit_functor_obstruction_audit.md"

class P3082Z12ContinuumLimitFunctorObstructionAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input_hash(self):
        self.assertEqual(self.payload["status"], "P3082_Z12_CONTINUUM_LIMIT_FUNCTOR_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P3081"])

    def test_finite_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 3)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["p3081_accepted_nonimported_dimension_action_unit_sources"], 0)
        self.assertEqual(cert["refinement_sizes"], 4)
        self.assertEqual(cert["mode_rows_per_refinement"], 5)
        self.assertEqual(cert["formal_refinement_spectral_rows"], 20)
        self.assertGreater(cert["formal_rows_below_error_0_1"], 0)
        self.assertLess(cert["max_final_refinement_error"], 0.25)
        self.assertEqual(cert["functor_candidates"], 5)
        self.assertEqual(cert["candidate_gate_rows"], 30)
        self.assertEqual(cert["accepted_nonimported_continuum_limit_functors"], 0)
        self.assertEqual(cert["satisfied_proof_obligations"], 4)

    def test_refinement_rows_and_candidate_aggregates(self):
        rows = self.payload["constructed_theoretical_objects"]["formal_refinement_spectral_rows"]
        self.assertEqual(sorted({row["cycle_size"] for row in rows}), [12, 24, 48, 96])
        self.assertEqual(sorted({row["mode"] for row in rows}), [1, 2, 3, 4, 5])
        mode1_errors = [row["absolute_error_after_imported_scaling"] for row in rows if row["mode"] == 1]
        self.assertLess(mode1_errors[-1], mode1_errors[0])
        aggs = self.payload["constructed_theoretical_objects"]["candidate_aggregate_certificate"]
        self.assertTrue(all(not row["accepted_nonimported_continuum_limit_functor"] for row in aggs))
        by_id = {row["candidate"]: row for row in aggs}
        self.assertEqual(by_id["fixed_z12_identity_functor"]["passed_gates"], 2)
        self.assertEqual(by_id["imported_standard_physics_manifold_template"]["passed_gates"], 5)

    def test_flags_recommendation_and_docs(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertTrue(self.payload["decision"]["positive_scoped_flags"]["formal_imported_fourier_convergence_witness_computed"])
        self.assertIn("Lorentz-signature", self.payload["decision"]["next_honest_step"])
        self.assertIn("P3082/S2032", MD.read_text(encoding="utf-8"))
        self.assertIn("P3082/S2032", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3082/S2032", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3082", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
