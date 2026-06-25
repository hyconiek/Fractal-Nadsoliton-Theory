import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3088_s2038_spectral_hamiltonian_time_evolution_obstruction_audit.py"
OUT = ROOT / "generated" / "p3088_s2038_spectral_hamiltonian_time_evolution_obstruction_audit.json"
MD = ROOT / "generated" / "p3088_s2038_spectral_hamiltonian_time_evolution_obstruction_audit.md"

class P3088SpectralHamiltonianTimeEvolutionObstructionAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input_hash(self):
        self.assertEqual(self.payload["status"], "P3088_SPECTRAL_HAMILTONIAN_TIME_EVOLUTION_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P3087"])

    def test_finite_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 3)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["p3087_accepted_nonimported_thermodynamic_ensemble_sources"], 0)
        self.assertEqual(cert["spectrum_rows"], 12)
        self.assertEqual(cert["spectral_hamiltonian_rows"], 12)
        self.assertEqual(cert["spectral_rows_with_energy_units"], 0)
        self.assertEqual(cert["formal_unitary_phase_rows"], 6)
        self.assertEqual(cert["unitarity_modulus_failures"], 0)
        self.assertEqual(cert["phase_rows_with_time_units"], 0)
        self.assertEqual(cert["phase_rows_with_action_units"], 0)
        self.assertEqual(cert["time_energy_scale_orbit_rows"], 18)
        self.assertEqual(cert["time_energy_scale_identity_failures"], 0)
        self.assertEqual(cert["canonical_time_or_energy_sources_exported"], 0)
        self.assertEqual(cert["hamiltonian_candidates"], 5)
        self.assertEqual(cert["candidate_gate_rows"], 30)
        self.assertEqual(cert["accepted_nonimported_hamiltonian_time_evolution_sources"], 0)
        self.assertEqual(cert["satisfied_proof_obligations"], 4)

    def test_hamiltonian_phase_and_scale_orbit(self):
        hrows = self.payload["constructed_theoretical_objects"]["spectral_hamiltonian_rows"]
        self.assertEqual(len(hrows), 12)
        self.assertTrue(all(row["finite_self_adjoint_real_symmetric_witness"] for row in hrows))
        self.assertAlmostEqual(hrows[0]["lambda_m"], 0.0)
        self.assertAlmostEqual(hrows[6]["lambda_m"], 4.0)
        phases = self.payload["constructed_theoretical_objects"]["formal_unitary_phase_rows"]
        self.assertTrue(all(row["formal_unitary_identity_holds"] for row in phases))
        self.assertTrue(all(not row["time_unit_attached"] for row in phases))
        orbit = self.payload["constructed_theoretical_objects"]["time_energy_scale_orbit_rows"]
        self.assertEqual(len(orbit), 18)
        self.assertTrue(all(row["scaled_energy_phase_matches_rescaled_time_phase"] for row in orbit))
        aggs = self.payload["constructed_theoretical_objects"]["candidate_aggregate_certificate"]
        self.assertTrue(all(not row["accepted_nonimported_hamiltonian_time_evolution_source"] for row in aggs))
        by_id = {row["candidate"]: row for row in aggs}
        self.assertEqual(by_id["z12_laplacian_spectral_multiplier"]["passed_gates"], 2)
        self.assertEqual(by_id["formal_schrodinger_lift_exp_minus_i_lambda_t"]["passed_gates"], 3)
        self.assertEqual(by_id["imported_quantum_mechanics_template"]["passed_gates"], 5)

    def test_flags_recommendation_and_docs(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertTrue(self.payload["decision"]["positive_scoped_flags"]["formal_unitary_phase_grid_computed"])
        self.assertIn("Born-rule probability-readout", self.payload["decision"]["next_honest_step"])
        self.assertIn("P3088/S2038", MD.read_text(encoding="utf-8"))
        self.assertIn("P3088/S2038", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3088/S2038", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3088", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
