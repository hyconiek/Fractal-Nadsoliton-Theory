import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3094_s2044_stress_energy_metric_response_obstruction_audit.py"
OUT = ROOT / "generated" / "p3094_s2044_stress_energy_metric_response_obstruction_audit.json"
MD = ROOT / "generated" / "p3094_s2044_stress_energy_metric_response_obstruction_audit.md"

class P3094StressEnergyMetricResponseObstructionAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input_hash(self):
        self.assertEqual(self.payload["status"], "P3094_STRESS_ENERGY_METRIC_RESPONSE_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P3093"])

    def test_finite_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 3)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["p3093_accepted_nonimported_ward_current_effective_charge_sources"], 0)
        self.assertEqual(cert["metric_variation_rows"], 144)
        self.assertEqual(cert["metric_variation_rows_with_metric_coupling"], 0)
        self.assertEqual(cert["metric_variation_rows_with_physical_stress_tensor"], 0)
        self.assertEqual(cert["graph_energy_quadratic_rows"], 12)
        self.assertEqual(cert["graph_energy_rows_with_action_measure"], 0)
        self.assertEqual(cert["graph_energy_rows_with_spacetime_metric"], 0)
        self.assertEqual(cert["spectral_pressure_like_rows"], 7)
        self.assertEqual(cert["spectral_pressure_rows_with_physical_pressure_unit"], 0)
        self.assertEqual(cert["spectral_pressure_rows_with_metric_response_law"], 0)
        self.assertEqual(cert["formal_stress_divergence_rows"], 12)
        self.assertEqual(cert["formal_stress_rows_with_covariant_conservation_law"], 0)
        self.assertEqual(cert["formal_stress_rows_with_empirical_field_response"], 0)
        self.assertEqual(cert["stress_candidates"], 5)
        self.assertEqual(cert["candidate_gate_rows"], 40)
        self.assertEqual(cert["accepted_nonimported_stress_energy_metric_response_sources"], 0)
        self.assertEqual(cert["satisfied_proof_obligations"], 4)

    def test_witnesses_and_candidates(self):
        objs = self.payload["constructed_theoretical_objects"]
        self.assertTrue(all(row["finite_metric_variation_witness"] for row in objs["metric_variation_rows"]))
        self.assertTrue(all(row["laplacian_quadratic_form_witness"] for row in objs["graph_energy_quadratic_rows"]))
        self.assertTrue(all(row["spectral_pressure_like_witness"] for row in objs["spectral_pressure_like_rows"]))
        self.assertTrue(all(row["formal_divergence_row_computed"] for row in objs["formal_stress_divergence_rows"]))
        aggs = objs["candidate_aggregate_certificate"]
        self.assertTrue(all(not row["accepted_nonimported_stress_energy_metric_response_source"] for row in aggs))
        by_id = {row["candidate"]: row for row in aggs}
        self.assertEqual(by_id["edge_weight_metric_variation_energy"]["passed_gates"], 3)
        self.assertEqual(by_id["imported_continuum_stress_tensor_template"]["passed_gates"], 6)
        self.assertEqual(by_id["imported_empirical_gravity_response_template"]["passed_gates"], 7)

    def test_flags_recommendation_and_docs(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertTrue(self.payload["decision"]["positive_scoped_flags"]["finite_metric_variation_rows_computed"])
        self.assertIn("dispersion/propagating-mode", self.payload["decision"]["next_honest_step"])
        self.assertIn("P3094/S2044", MD.read_text(encoding="utf-8"))
        self.assertIn("P3094/S2044", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3094/S2044", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3094", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
