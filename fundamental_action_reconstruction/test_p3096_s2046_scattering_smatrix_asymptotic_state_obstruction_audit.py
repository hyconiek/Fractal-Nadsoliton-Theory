import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3096_s2046_scattering_smatrix_asymptotic_state_obstruction_audit.py"
OUT = ROOT / "generated" / "p3096_s2046_scattering_smatrix_asymptotic_state_obstruction_audit.json"
MD = ROOT / "generated" / "p3096_s2046_scattering_smatrix_asymptotic_state_obstruction_audit.md"

class P3096ScatteringSMatrixAsymptoticStateObstructionAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input_hash(self):
        self.assertEqual(self.payload["status"], "P3096_SCATTERING_SMATRIX_ASYMPTOTIC_STATE_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P3095"])

    def test_finite_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 3)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["p3095_accepted_nonimported_dispersion_propagating_observable_sources"], 0)
        self.assertEqual(cert["finite_channel_rows"], 12)
        self.assertEqual(cert["channel_rows_with_in_state"], 0)
        self.assertEqual(cert["channel_rows_with_out_state"], 0)
        self.assertEqual(cert["channel_rows_with_asymptotic_region"], 0)
        self.assertEqual(cert["born_transition_amplitude_rows"], 144)
        self.assertEqual(cert["transition_rows_with_unit_normalized_amplitude"], 0)
        self.assertEqual(cert["transition_rows_with_cross_section_semantics"], 0)
        self.assertEqual(cert["s_matrix_unitarity_proxy_rows"], 3)
        self.assertEqual(cert["s_proxy_rows_with_exact_unitarity"], 0)
        self.assertEqual(cert["s_proxy_rows_with_physical_scattering_operator"], 0)
        self.assertEqual(cert["cross_section_proxy_rows"], 12)
        self.assertEqual(cert["cross_section_rows_with_area_unit"], 0)
        self.assertEqual(cert["cross_section_rows_with_detector_semantics"], 0)
        self.assertEqual(cert["scattering_candidates"], 5)
        self.assertEqual(cert["candidate_gate_rows"], 40)
        self.assertEqual(cert["accepted_nonimported_scattering_smatrix_sources"], 0)
        self.assertEqual(cert["satisfied_proof_obligations"], 4)

    def test_witnesses_and_candidates(self):
        objs = self.payload["constructed_theoretical_objects"]
        self.assertTrue(all(row["finite_channel_basis_witness"] for row in objs["finite_channel_rows"]))
        self.assertTrue(all(row["finite_transition_amplitude_witness"] for row in objs["born_transition_amplitude_rows"]))
        self.assertTrue(all(row["finite_s_matrix_proxy_witness"] for row in objs["s_matrix_unitarity_proxy_rows"]))
        self.assertTrue(all(row["finite_cross_section_proxy_witness"] for row in objs["cross_section_proxy_rows"]))
        aggs = objs["candidate_aggregate_certificate"]
        self.assertTrue(all(not row["accepted_nonimported_scattering_smatrix_source"] for row in aggs))
        by_id = {row["candidate"]: row for row in aggs}
        self.assertEqual(by_id["finite_fourier_channel_catalog"]["passed_gates"], 1)
        self.assertEqual(by_id["finite_born_transition_matrix"]["passed_gates"], 2)
        self.assertEqual(by_id["imported_empirical_detector_scattering_template"]["passed_gates"], 7)

    def test_flags_recommendation_and_docs(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertTrue(self.payload["decision"]["positive_scoped_flags"]["born_transition_amplitude_rows_computed"])
        self.assertIn("thermodynamic-radiation/blackbody-spectrum", self.payload["decision"]["next_honest_step"])
        self.assertIn("P3096/S2046", MD.read_text(encoding="utf-8"))
        self.assertIn("P3096/S2046", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3096/S2046", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3096", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
