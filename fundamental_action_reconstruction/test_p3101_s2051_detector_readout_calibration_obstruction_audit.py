import json, subprocess, sys, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3101_s2051_detector_readout_calibration_obstruction_audit.py"
OUT = ROOT / "generated" / "p3101_s2051_detector_readout_calibration_obstruction_audit.json"
MD = ROOT / "generated" / "p3101_s2051_detector_readout_calibration_obstruction_audit.md"

class P3101DetectorReadoutCalibrationObstructionAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_predecessor_hash(self):
        self.assertEqual(self.payload["status"], "P3101_DETECTOR_READOUT_CALIBRATION_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertTrue(self.payload["input_hashes"]["P3100"])

    def test_certificate_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["p3100_accepted_nonimported_bath_preparation_sources"], 0)
        self.assertEqual(cert["detector_response_rows"], 48)
        self.assertEqual(cert["detector_response_rows_with_empirical_detector"], 0)
        self.assertEqual(cert["calibration_orbit_rows"], 60)
        self.assertEqual(cert["calibration_rows_with_canonical_unit"], 0)
        self.assertEqual(cert["threshold_classifier_rows"], 48)
        self.assertEqual(cert["threshold_rows_with_physical_threshold_source"], 0)
        self.assertEqual(cert["noise_stability_rows"], 144)
        self.assertEqual(cert["noise_rows_with_physical_noise_model"], 0)
        self.assertEqual(cert["detector_candidates"], 5)
        self.assertEqual(cert["required_gates"], 8)
        self.assertEqual(cert["candidate_gate_rows"], 40)
        self.assertEqual(cert["accepted_nonimported_detector_sources"], 0)

    def test_witness_rows_are_formal_and_unsourced(self):
        objs = self.payload["constructed_theoretical_objects"]
        self.assertTrue(all(row["finite_response_map_witness"] for row in objs["detector_response_rows"]))
        self.assertTrue(all(not row["empirical_detector_attached"] for row in objs["detector_response_rows"]))
        self.assertTrue(all(row["scale_orbit_witness"] for row in objs["calibration_orbit_rows"]))
        self.assertTrue(all(not row["canonical_unit_fixed"] for row in objs["calibration_orbit_rows"]))
        self.assertTrue(all(row["threshold_classifier_witness"] for row in objs["threshold_classifier_rows"]))
        self.assertTrue(all(not row["physical_threshold_source_attached"] for row in objs["threshold_classifier_rows"]))
        self.assertTrue(all(row["noise_stability_witness"] for row in objs["noise_stability_rows"]))
        self.assertTrue(all(not row["physical_noise_model_attached"] for row in objs["noise_stability_rows"]))

    def test_negative_flags_and_recommendation(self):
        decision = self.payload["decision"]
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))
        self.assertIn("Born-rule/probability-measure", decision["next_honest_step"])
        self.assertIn("P3101/S2051", MD.read_text(encoding="utf-8"))
        self.assertIn("P3101/S2051", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3101/S2051", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3101", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
