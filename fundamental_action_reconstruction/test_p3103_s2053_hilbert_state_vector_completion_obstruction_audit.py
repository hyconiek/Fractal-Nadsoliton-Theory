import json, subprocess, sys, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3103_s2053_hilbert_state_vector_completion_obstruction_audit.py"
OUT = ROOT / "generated" / "p3103_s2053_hilbert_state_vector_completion_obstruction_audit.json"
MD = ROOT / "generated" / "p3103_s2053_hilbert_state_vector_completion_obstruction_audit.md"

class P3103HilbertStateVectorCompletionObstructionAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_predecessor_hash(self):
        self.assertEqual(self.payload["status"], "P3103_HILBERT_STATE_VECTOR_COMPLETION_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertTrue(self.payload["input_hashes"]["P3102"])

    def test_certificate_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["p3102_accepted_nonimported_probability_sources"], 0)
        self.assertEqual(cert["inner_product_rows"], 16)
        self.assertEqual(cert["inner_product_rows_with_physical_units"], 0)
        self.assertEqual(cert["observable_expectation_rows"], 12)
        self.assertEqual(cert["observable_rows_with_physical_units"], 0)
        self.assertEqual(cert["unitary_orbit_rows"], 4)
        self.assertEqual(cert["unitary_rows_with_physical_time"], 0)
        self.assertEqual(cert["preparation_readout_rows"], 48)
        self.assertEqual(cert["preparation_rows_with_physical_source"], 0)
        self.assertEqual(cert["hilbert_candidates"], 5)
        self.assertEqual(cert["required_gates"], 8)
        self.assertEqual(cert["candidate_gate_rows"], 40)
        self.assertEqual(cert["accepted_nonimported_hilbert_sources"], 0)

    def test_witness_rows_are_formal_and_unsourced(self):
        objs = self.payload["constructed_theoretical_objects"]
        self.assertTrue(all(row["hermitian_partner_expected"] for row in objs["inner_product_rows"]))
        self.assertTrue(all(not row["physical_unit_attached"] for row in objs["inner_product_rows"]))
        self.assertTrue(all(row["self_adjoint_proxy_witness"] for row in objs["observable_expectation_rows"]))
        self.assertTrue(all(not row["observable_physical_units_attached"] for row in objs["observable_expectation_rows"]))
        self.assertTrue(all(row["unitarity_proxy_witness"] for row in objs["unitary_orbit_rows"]))
        self.assertTrue(all(not row["physical_time_parameter_attached"] for row in objs["unitary_orbit_rows"]))
        self.assertTrue(all(row["preparation_map_proxy"] for row in objs["preparation_readout_rows"]))
        self.assertTrue(all(not row["physical_preparation_source_attached"] for row in objs["preparation_readout_rows"]))

    def test_negative_flags_and_recommendation(self):
        decision = self.payload["decision"]
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))
        self.assertIn("spectral-triple/geometry-interface", decision["next_honest_step"])
        self.assertIn("P3103/S2053", MD.read_text(encoding="utf-8"))
        self.assertIn("P3103/S2053", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3103/S2053", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3103", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
