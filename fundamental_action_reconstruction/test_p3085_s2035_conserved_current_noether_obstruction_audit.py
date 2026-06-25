import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3085_s2035_conserved_current_noether_obstruction_audit.py"
OUT = ROOT / "generated" / "p3085_s2035_conserved_current_noether_obstruction_audit.json"
MD = ROOT / "generated" / "p3085_s2035_conserved_current_noether_obstruction_audit.md"

class P3085ConservedCurrentNoetherObstructionAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input_hash(self):
        self.assertEqual(self.payload["status"], "P3085_CONSERVED_CURRENT_NOETHER_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P3084"])

    def test_finite_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 3)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["p3084_accepted_nonimported_gauge_representation_sources"], 0)
        self.assertEqual(cert["fourier_link_current_rows"], 12)
        self.assertEqual(cert["formal_divergence_free_fourier_rows"], 12)
        self.assertEqual(cert["fourier_rows_with_unit_bearing_current"], 0)
        self.assertEqual(cert["finite_continuity_matrix_rows"], 144)
        self.assertEqual(cert["continuity_matrix_failures"], 0)
        self.assertEqual(cert["binary_real_profile_control_rows"], 4096)
        self.assertEqual(cert["binary_rows_with_phase_degree"], 0)
        self.assertEqual(cert["binary_rows_with_charge_density"], 0)
        self.assertEqual(cert["current_candidates"], 5)
        self.assertEqual(cert["candidate_gate_rows"], 30)
        self.assertEqual(cert["accepted_nonimported_conserved_current_sources"], 0)
        self.assertEqual(cert["satisfied_proof_obligations"], 4)

    def test_current_rows_and_candidate_aggregates(self):
        rows = {row["mode_m"]: row for row in self.payload["constructed_theoretical_objects"]["fourier_link_current_rows"]}
        self.assertAlmostEqual(rows[0]["link_current_im_conj_psi_i_psi_ip1"], 0.0)
        self.assertAlmostEqual(rows[3]["link_current_im_conj_psi_i_psi_ip1"], 1.0)
        self.assertTrue(all(row["formal_continuity_witness"] for row in rows.values()))
        aggs = self.payload["constructed_theoretical_objects"]["candidate_aggregate_certificate"]
        self.assertTrue(all(not row["accepted_nonimported_conserved_current_source"] for row in aggs))
        by_id = {row["candidate"]: row for row in aggs}
        self.assertEqual(by_id["real_dirichlet_energy_density"]["passed_gates"], 1)
        self.assertEqual(by_id["z12_fourier_link_current_witness"]["passed_gates"], 2)
        self.assertEqual(by_id["imported_continuum_noether_template"]["passed_gates"], 5)

    def test_flags_recommendation_and_docs(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertTrue(self.payload["decision"]["positive_scoped_flags"]["finite_continuity_rows_verified_formally"])
        self.assertIn("empirical-readout/observable-calibration", self.payload["decision"]["next_honest_step"])
        self.assertIn("P3085/S2035", MD.read_text(encoding="utf-8"))
        self.assertIn("P3085/S2035", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3085/S2035", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3085", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
