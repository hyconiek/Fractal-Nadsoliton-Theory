import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3087_s2037_thermodynamic_statistical_ensemble_obstruction_audit.py"
OUT = ROOT / "generated" / "p3087_s2037_thermodynamic_statistical_ensemble_obstruction_audit.json"
MD = ROOT / "generated" / "p3087_s2037_thermodynamic_statistical_ensemble_obstruction_audit.md"

class P3087ThermodynamicStatisticalEnsembleObstructionAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input_hash(self):
        self.assertEqual(self.payload["status"], "P3087_THERMODYNAMIC_STATISTICAL_ENSEMBLE_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P3086"])

    def test_finite_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 3)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["p3086_accepted_nonimported_empirical_observable_sources"], 0)
        self.assertEqual(cert["spectrum_rows"], 12)
        self.assertEqual(cert["microcanonical_degeneracy_rows"], 7)
        self.assertEqual(cert["formal_partition_function_rows"], 7)
        self.assertEqual(cert["partition_rows_with_temperature_units"], 0)
        self.assertEqual(cert["partition_rows_with_energy_units"], 0)
        self.assertEqual(cert["energy_scale_beta_orbit_rows"], 21)
        self.assertEqual(cert["scale_beta_identity_failures"], 0)
        self.assertEqual(cert["canonical_temperature_sources_exported"], 0)
        self.assertEqual(cert["ensemble_candidates"], 5)
        self.assertEqual(cert["candidate_gate_rows"], 30)
        self.assertEqual(cert["accepted_nonimported_thermodynamic_ensemble_sources"], 0)
        self.assertEqual(cert["satisfied_proof_obligations"], 4)

    def test_partition_degeneracy_and_scale_orbit(self):
        degeneracies = self.payload["constructed_theoretical_objects"]["microcanonical_degeneracy_rows"]
        self.assertEqual(sum(row["degeneracy"] for row in degeneracies), 12)
        self.assertEqual(degeneracies[0]["energy_label"], 0.0)
        self.assertEqual(degeneracies[-1]["energy_label"], 4.0)
        partitions = {row["dimensionless_beta"]: row for row in self.payload["constructed_theoretical_objects"]["formal_partition_function_rows"]}
        self.assertAlmostEqual(partitions[0.0]["partition_function_Z"], 12.0)
        self.assertTrue(partitions[4.0]["partition_function_Z"] < partitions[0.0]["partition_function_Z"])
        orbit = self.payload["constructed_theoretical_objects"]["energy_scale_beta_orbit_rows"]
        self.assertTrue(all(row["scale_beta_compensation_identity_holds"] for row in orbit))
        aggs = self.payload["constructed_theoretical_objects"]["candidate_aggregate_certificate"]
        self.assertTrue(all(not row["accepted_nonimported_thermodynamic_ensemble_source"] for row in aggs))
        by_id = {row["candidate"]: row for row in aggs}
        self.assertEqual(by_id["finite_z12_laplacian_boltzmann_weights"]["passed_gates"], 2)
        self.assertEqual(by_id["imported_boltzmann_gibbs_template"]["passed_gates"], 5)

    def test_flags_recommendation_and_docs(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertTrue(self.payload["decision"]["positive_scoped_flags"]["formal_partition_grid_computed"])
        self.assertIn("spectral-to-Hamiltonian/time-evolution", self.payload["decision"]["next_honest_step"])
        self.assertIn("P3087/S2037", MD.read_text(encoding="utf-8"))
        self.assertIn("P3087/S2037", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3087/S2037", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3087", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
