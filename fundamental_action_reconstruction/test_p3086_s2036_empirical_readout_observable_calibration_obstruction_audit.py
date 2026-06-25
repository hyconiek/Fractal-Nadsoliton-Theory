import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3086_s2036_empirical_readout_observable_calibration_obstruction_audit.py"
OUT = ROOT / "generated" / "p3086_s2036_empirical_readout_observable_calibration_obstruction_audit.json"
MD = ROOT / "generated" / "p3086_s2036_empirical_readout_observable_calibration_obstruction_audit.md"

class P3086EmpiricalReadoutObservableCalibrationObstructionAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input_hash(self):
        self.assertEqual(self.payload["status"], "P3086_EMPIRICAL_READOUT_OBSERVABLE_CALIBRATION_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P3085"])

    def test_finite_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 3)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["p3085_accepted_nonimported_conserved_current_sources"], 0)
        self.assertEqual(cert["spectrum_rows"], 12)
        self.assertEqual(cert["spectrum_rows_with_unit_calibration"], 0)
        self.assertEqual(cert["distinct_spectral_gap_rows"], 7)
        self.assertEqual(cert["scale_orbit_control_rows"], 5)
        self.assertEqual(cert["canonical_scale_sources_exported"], 0)
        self.assertEqual(cert["observable_candidates"], 5)
        self.assertEqual(cert["candidate_gate_rows"], 30)
        self.assertEqual(cert["calibration_attempt_rows"], 30)
        self.assertEqual(cert["accepted_nonimported_empirical_observable_sources"], 0)
        self.assertEqual(cert["satisfied_proof_obligations"], 4)

    def test_spectrum_scale_controls_and_candidate_aggregates(self):
        rows = {row["mode_m"]: row for row in self.payload["constructed_theoretical_objects"]["z12_laplacian_spectrum_rows"]}
        self.assertAlmostEqual(rows[0]["laplacian_eigenvalue"], 0.0)
        self.assertAlmostEqual(rows[6]["laplacian_eigenvalue"], 4.0)
        self.assertTrue(all(not row["unit_calibrated"] for row in rows.values()))
        scales = self.payload["constructed_theoretical_objects"]["scale_orbit_control_rows"]
        self.assertEqual([row["scale_factor"] for row in scales], [0.25, 0.5, 1.0, 2.0, 4.0])
        self.assertTrue(all(row["ratios_preserved"] for row in scales))
        aggs = self.payload["constructed_theoretical_objects"]["candidate_aggregate_certificate"]
        self.assertTrue(all(not row["accepted_nonimported_empirical_observable_source"] for row in aggs))
        by_id = {row["candidate"]: row for row in aggs}
        self.assertEqual(by_id["z12_dirichlet_scalar_energy_proxy"]["passed_gates"], 2)
        self.assertEqual(by_id["z12_laplacian_spectrum_proxy"]["passed_gates"], 2)
        self.assertEqual(by_id["imported_dimensionful_meter_calibration_template"]["passed_gates"], 5)

    def test_flags_recommendation_and_docs(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertTrue(self.payload["decision"]["positive_scoped_flags"]["scale_orbit_controls_executed"])
        self.assertIn("thermodynamic/statistical-ensemble", self.payload["decision"]["next_honest_step"])
        self.assertIn("P3086/S2036", MD.read_text(encoding="utf-8"))
        self.assertIn("P3086/S2036", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3086/S2036", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3086", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
