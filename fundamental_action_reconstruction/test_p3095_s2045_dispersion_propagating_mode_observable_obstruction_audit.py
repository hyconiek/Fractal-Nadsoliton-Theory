import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3095_s2045_dispersion_propagating_mode_observable_obstruction_audit.py"
OUT = ROOT / "generated" / "p3095_s2045_dispersion_propagating_mode_observable_obstruction_audit.json"
MD = ROOT / "generated" / "p3095_s2045_dispersion_propagating_mode_observable_obstruction_audit.md"

class P3095DispersionPropagatingModeObservableObstructionAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input_hash(self):
        self.assertEqual(self.payload["status"], "P3095_DISPERSION_PROPAGATING_MODE_OBSERVABLE_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P3094"])

    def test_finite_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 3)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["p3094_accepted_nonimported_stress_energy_metric_response_sources"], 0)
        self.assertEqual(cert["dispersion_group_velocity_rows"], 12)
        self.assertEqual(cert["dispersion_rows_with_spacetime_speed"], 0)
        self.assertEqual(cert["dispersion_rows_with_observed_light_semantics"], 0)
        self.assertEqual(cert["mode_packet_evolution_rows"], 12)
        self.assertEqual(cert["packet_rows_with_physical_time_unit"], 0)
        self.assertEqual(cert["packet_rows_with_detector_observable"], 0)
        self.assertEqual(cert["green_pole_catalog_rows"], 48)
        self.assertEqual(cert["green_pole_rows_with_asymptotic_state"], 0)
        self.assertEqual(cert["green_pole_rows_with_scattering_readout"], 0)
        self.assertEqual(cert["energy_flux_proxy_rows"], 12)
        self.assertEqual(cert["energy_flux_rows_with_physical_flux_unit"], 0)
        self.assertEqual(cert["energy_flux_rows_with_radiation_interface"], 0)
        self.assertEqual(cert["propagation_candidates"], 5)
        self.assertEqual(cert["candidate_gate_rows"], 40)
        self.assertEqual(cert["accepted_nonimported_dispersion_propagating_observable_sources"], 0)
        self.assertEqual(cert["satisfied_proof_obligations"], 4)

    def test_witnesses_and_candidates(self):
        objs = self.payload["constructed_theoretical_objects"]
        self.assertTrue(all(row["finite_dispersion_witness"] for row in objs["dispersion_group_velocity_rows"]))
        self.assertTrue(all(row["formal_mode_packet_evolution_witness"] for row in objs["mode_packet_evolution_rows"]))
        self.assertTrue(all(row["green_pole_catalog_witness"] for row in objs["green_pole_catalog_rows"]))
        self.assertTrue(all(row["stress_energy_link_witness"] for row in objs["energy_flux_proxy_rows"]))
        aggs = objs["candidate_aggregate_certificate"]
        self.assertTrue(all(not row["accepted_nonimported_dispersion_propagating_observable_source"] for row in aggs))
        by_id = {row["candidate"]: row for row in aggs}
        self.assertEqual(by_id["finite_laplacian_dispersion_curve"]["passed_gates"], 2)
        self.assertEqual(by_id["formal_unitary_mode_packet_evolution"]["passed_gates"], 3)
        self.assertEqual(by_id["imported_observed_light_radiation_template"]["passed_gates"], 7)

    def test_flags_recommendation_and_docs(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertTrue(self.payload["decision"]["positive_scoped_flags"]["finite_dispersion_group_velocity_rows_computed"])
        self.assertIn("scattering/S-matrix", self.payload["decision"]["next_honest_step"])
        self.assertIn("P3095/S2045", MD.read_text(encoding="utf-8"))
        self.assertIn("P3095/S2045", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3095/S2045", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3095", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
