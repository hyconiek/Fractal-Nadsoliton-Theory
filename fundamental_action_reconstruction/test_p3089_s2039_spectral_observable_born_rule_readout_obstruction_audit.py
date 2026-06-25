import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3089_s2039_spectral_observable_born_rule_readout_obstruction_audit.py"
OUT = ROOT / "generated" / "p3089_s2039_spectral_observable_born_rule_readout_obstruction_audit.json"
MD = ROOT / "generated" / "p3089_s2039_spectral_observable_born_rule_readout_obstruction_audit.md"

class P3089SpectralObservableBornRuleReadoutObstructionAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input_hash(self):
        self.assertEqual(self.payload["status"], "P3089_SPECTRAL_OBSERVABLE_BORN_RULE_READOUT_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P3088"])

    def test_finite_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 3)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["p3088_accepted_nonimported_hamiltonian_time_evolution_sources"], 0)
        self.assertEqual(cert["nonzero_binary_profile_probability_census_rows"], 4095)
        self.assertEqual(cert["normalized_probability_failures"], 0)
        self.assertEqual(cert["negative_probability_failures"], 0)
        self.assertEqual(cert["rows_with_born_rule_source_attached"], 0)
        self.assertEqual(cert["rows_with_empirical_frequency_readout_attached"], 0)
        self.assertEqual(cert["translation_orbit_probability_rows"], 12)
        self.assertEqual(cert["translation_invariance_failures"], 0)
        self.assertEqual(cert["source_origin_localizer_rows"], 0)
        self.assertEqual(cert["formal_time_probability_conservation_rows"], 6)
        self.assertEqual(cert["formal_time_probability_conservation_failures"], 0)
        self.assertEqual(cert["time_rows_with_measurement_readout"], 0)
        self.assertEqual(cert["born_rule_candidates"], 5)
        self.assertEqual(cert["candidate_gate_rows"], 30)
        self.assertEqual(cert["accepted_nonimported_born_rule_probability_readout_sources"], 0)
        self.assertEqual(cert["satisfied_proof_obligations"], 4)

    def test_probability_census_and_controls(self):
        census = self.payload["constructed_theoretical_objects"]["nonzero_binary_profile_probability_census_rows"]
        self.assertEqual(census[0]["mask"], 1)
        self.assertEqual(census[-1]["mask"], 4095)
        self.assertTrue(all(row["normalized_probability_witness"] for row in census))
        self.assertTrue(all(row["nonnegative_probability_witness"] for row in census))
        shifts = self.payload["constructed_theoretical_objects"]["source_translation_orbit_probability_rows"]
        self.assertTrue(all(row["power_spectrum_translation_invariant"] for row in shifts))
        self.assertTrue(all(not row["source_origin_localized_by_probability"] for row in shifts))
        time_rows = self.payload["constructed_theoretical_objects"]["formal_time_probability_conservation_rows"]
        self.assertTrue(all(row["formal_probability_conservation_holds"] for row in time_rows))
        aggs = self.payload["constructed_theoretical_objects"]["candidate_aggregate_certificate"]
        self.assertTrue(all(not row["accepted_nonimported_born_rule_probability_readout_source"] for row in aggs))
        by_id = {row["candidate"]: row for row in aggs}
        self.assertEqual(by_id["z12_normalized_fourier_power_spectrum"]["passed_gates"], 2)
        self.assertEqual(by_id["p3088_formal_unitary_probability_conservation"]["passed_gates"], 2)
        self.assertEqual(by_id["imported_empirical_detector_frequency_template"]["passed_gates"], 5)

    def test_flags_recommendation_and_docs(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertTrue(self.payload["decision"]["positive_scoped_flags"]["formal_time_probability_conservation_verified"])
        self.assertIn("Green-function response", self.payload["decision"]["next_honest_step"])
        self.assertIn("P3089/S2039", MD.read_text(encoding="utf-8"))
        self.assertIn("P3089/S2039", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3089/S2039", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3089", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
