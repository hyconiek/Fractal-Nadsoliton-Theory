import json, subprocess, sys, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3102_s2052_born_probability_measure_empirical_readout_obstruction_audit.py"
OUT = ROOT / "generated" / "p3102_s2052_born_probability_measure_empirical_readout_obstruction_audit.json"
MD = ROOT / "generated" / "p3102_s2052_born_probability_measure_empirical_readout_obstruction_audit.md"

class P3102BornProbabilityMeasureEmpiricalReadoutObstructionAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_predecessor_hash(self):
        self.assertEqual(self.payload["status"], "P3102_BORN_PROBABILITY_MEASURE_EMPIRICAL_READOUT_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertTrue(self.payload["input_hashes"]["P3101"])

    def test_certificate_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["p3101_accepted_nonimported_detector_sources"], 0)
        self.assertEqual(cert["probability_measure_rows"], 48)
        self.assertEqual(cert["probability_rows_with_born_rule_source"], 0)
        self.assertEqual(cert["event_additivity_rows"], 16)
        self.assertEqual(cert["event_rows_with_physical_semantics"], 0)
        self.assertEqual(cert["basis_orbit_rows"], 72)
        self.assertEqual(cert["basis_rows_with_canonical_measurement_basis"], 0)
        self.assertEqual(cert["frequency_proxy_rows"], 16)
        self.assertEqual(cert["frequency_rows_with_empirical_readout"], 0)
        self.assertEqual(cert["probability_candidates"], 5)
        self.assertEqual(cert["required_gates"], 8)
        self.assertEqual(cert["candidate_gate_rows"], 40)
        self.assertEqual(cert["accepted_nonimported_probability_sources"], 0)

    def test_witness_rows_are_formal_and_unsourced(self):
        objs = self.payload["constructed_theoretical_objects"]
        self.assertTrue(all(row["normalized_measure_witness"] for row in objs["probability_measure_rows"]))
        self.assertTrue(all(not row["born_rule_source_attached"] for row in objs["probability_measure_rows"]))
        self.assertTrue(all(row["finite_additivity_witness"] for row in objs["event_additivity_rows"]))
        self.assertTrue(all(not row["physical_event_semantics_attached"] for row in objs["event_additivity_rows"]))
        self.assertTrue(all(row["basis_orbit_witness"] for row in objs["basis_orbit_rows"]))
        self.assertTrue(all(not row["canonical_measurement_basis_fixed"] for row in objs["basis_orbit_rows"]))
        self.assertTrue(all(row["frequency_proxy_witness"] for row in objs["frequency_proxy_rows"]))
        self.assertTrue(all(not row["empirical_readout_attached"] for row in objs["frequency_proxy_rows"]))

    def test_negative_flags_and_recommendation(self):
        decision = self.payload["decision"]
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))
        self.assertIn("Hilbert-space/state-vector", decision["next_honest_step"])
        self.assertIn("P3102/S2052", MD.read_text(encoding="utf-8"))
        self.assertIn("P3102/S2052", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3102/S2052", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3102", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
