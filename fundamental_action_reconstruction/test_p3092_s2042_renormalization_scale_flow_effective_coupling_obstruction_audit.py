import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3092_s2042_renormalization_scale_flow_effective_coupling_obstruction_audit.py"
OUT = ROOT / "generated" / "p3092_s2042_renormalization_scale_flow_effective_coupling_obstruction_audit.json"
MD = ROOT / "generated" / "p3092_s2042_renormalization_scale_flow_effective_coupling_obstruction_audit.md"

class P3092RenormalizationScaleFlowObstructionAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input_hash(self):
        self.assertEqual(self.payload["status"], "P3092_RENORMALIZATION_SCALE_FLOW_EFFECTIVE_COUPLING_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P3091"])

    def test_finite_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 3)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["p3091_accepted_nonimported_spectral_action_effective_action_sources"], 0)
        self.assertEqual(cert["logdet_scale_dependence_rows"], 7)
        self.assertEqual(cert["logdet_rows_with_physical_rg_scale"], 0)
        self.assertEqual(cert["logdet_rows_with_action_unit"], 0)
        self.assertEqual(cert["green_beta_like_rows"], 6)
        self.assertEqual(cert["green_beta_like_rows_with_sourced_beta_function"], 0)
        self.assertEqual(cert["green_beta_like_rows_with_empirical_matching"], 0)
        self.assertEqual(cert["coupling_rescaling_orbit_rows"], 28)
        self.assertEqual(cert["coupling_rows_with_unit_normalized_running_coupling"], 0)
        self.assertEqual(cert["coupling_rows_with_physical_normalization"], 0)
        self.assertEqual(cert["rg_candidates"], 5)
        self.assertEqual(cert["candidate_gate_rows"], 35)
        self.assertEqual(cert["accepted_nonimported_renormalization_scale_flow_sources"], 0)
        self.assertEqual(cert["satisfied_proof_obligations"], 4)

    def test_witnesses_and_candidates(self):
        objs = self.payload["constructed_theoretical_objects"]
        self.assertTrue(all(row["finite_scale_dependence_witness"] for row in objs["logdet_scale_dependence_rows"]))
        self.assertTrue(all(row["beta_like_witness_computed"] for row in objs["green_effective_coupling_beta_like_rows"]))
        self.assertTrue(all(row["coupling_rescaling_witness"] for row in objs["coupling_rescaling_orbit_rows"]))
        aggs = objs["candidate_aggregate_certificate"]
        self.assertTrue(all(not row["accepted_nonimported_renormalization_scale_flow_source"] for row in aggs))
        by_id = {row["candidate"]: row for row in aggs}
        self.assertEqual(by_id["finite_logdet_mass_scale_derivative"]["passed_gates"], 2)
        self.assertEqual(by_id["green_ratio_effective_coupling_orbit"]["passed_gates"], 1)
        self.assertEqual(by_id["imported_empirical_running_coupling_template"]["passed_gates"], 6)

    def test_flags_recommendation_and_docs(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertTrue(self.payload["decision"]["positive_scoped_flags"]["green_beta_like_finite_differences_computed"])
        self.assertIn("Ward-identity/symmetry-current", self.payload["decision"]["next_honest_step"])
        self.assertIn("P3092/S2042", MD.read_text(encoding="utf-8"))
        self.assertIn("P3092/S2042", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3092/S2042", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3092", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
