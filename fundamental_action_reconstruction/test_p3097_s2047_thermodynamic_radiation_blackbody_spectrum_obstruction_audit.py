import json, subprocess, sys, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3097_s2047_thermodynamic_radiation_blackbody_spectrum_obstruction_audit.py"
OUT = ROOT / "generated" / "p3097_s2047_thermodynamic_radiation_blackbody_spectrum_obstruction_audit.json"
MD = ROOT / "generated" / "p3097_s2047_thermodynamic_radiation_blackbody_spectrum_obstruction_audit.md"

class P3097ThermodynamicRadiationBlackbodySpectrumObstructionAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_predecessor_hash(self):
        self.assertEqual(self.payload["status"], "P3097_THERMODYNAMIC_RADIATION_BLACKBODY_SPECTRUM_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertTrue(self.payload["input_hashes"]["P3096"])

    def test_certificate_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["p3096_accepted_nonimported_scattering_smatrix_sources"], 0)
        self.assertEqual(cert["finite_mode_count_rows"], 7)
        self.assertEqual(cert["formal_partition_weight_rows"], 4)
        self.assertEqual(cert["radiation_spectrum_proxy_rows"], 28)
        self.assertEqual(cert["temperature_energy_scale_orbit_rows"], 3)
        self.assertEqual(cert["radiation_candidates"], 5)
        self.assertEqual(cert["required_gates"], 7)
        self.assertEqual(cert["candidate_gate_rows"], 35)
        self.assertEqual(cert["accepted_nonimported_thermodynamic_radiation_sources"], 0)

    def test_witness_rows_are_bounded_and_formal(self):
        objs = self.payload["constructed_theoretical_objects"]
        self.assertEqual(len(objs["finite_mode_count_rows"]), 7)
        self.assertTrue(all(row["finite_mode_count_witness"] for row in objs["finite_mode_count_rows"]))
        self.assertTrue(all(not row["temperature_unit_attached"] for row in objs["finite_mode_count_rows"]))
        self.assertTrue(all(row["dimensionless_partition_witness"] for row in objs["formal_partition_weight_rows"]))
        self.assertTrue(all(row["radiation_spectrum_proxy_witness"] for row in objs["radiation_spectrum_proxy_rows"]))
        self.assertTrue(all(not row["observed_light_semantics_attached"] for row in objs["radiation_spectrum_proxy_rows"]))
        self.assertTrue(all(row["scale_orbit_witness"] for row in objs["temperature_energy_scale_orbit_rows"]))

    def test_negative_flags_and_recommendation(self):
        decision = self.payload["decision"]
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))
        self.assertIn("KMS/detailed-balance thermal-state", decision["next_honest_step"])
        self.assertIn("P3097/S2047", MD.read_text(encoding="utf-8"))
        self.assertIn("P3097/S2047", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3097/S2047", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3097", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
