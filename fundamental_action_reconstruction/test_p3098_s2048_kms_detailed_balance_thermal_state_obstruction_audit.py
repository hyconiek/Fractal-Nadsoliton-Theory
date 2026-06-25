import json, subprocess, sys, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3098_s2048_kms_detailed_balance_thermal_state_obstruction_audit.py"
OUT = ROOT / "generated" / "p3098_s2048_kms_detailed_balance_thermal_state_obstruction_audit.json"
MD = ROOT / "generated" / "p3098_s2048_kms_detailed_balance_thermal_state_obstruction_audit.md"

class P3098KMSDetailedBalanceThermalStateObstructionAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_predecessor_hash(self):
        self.assertEqual(self.payload["status"], "P3098_KMS_DETAILED_BALANCE_THERMAL_STATE_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertTrue(self.payload["input_hashes"]["P3097"])

    def test_certificate_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["p3097_accepted_nonimported_thermodynamic_radiation_sources"], 0)
        self.assertEqual(cert["gibbs_weight_rows"], 48)
        self.assertEqual(cert["detailed_balance_transition_rows"], 528)
        self.assertEqual(cert["detailed_balance_rows_with_zero_residual"], 528)
        self.assertEqual(cert["kms_periodicity_proxy_rows"], 16)
        self.assertEqual(cert["fluctuation_dissipation_proxy_rows"], 40)
        self.assertEqual(cert["thermal_state_candidates"], 5)
        self.assertEqual(cert["required_gates"], 8)
        self.assertEqual(cert["candidate_gate_rows"], 40)
        self.assertEqual(cert["accepted_nonimported_kms_thermal_state_sources"], 0)

    def test_witness_rows_are_exact_but_unsourced(self):
        objs = self.payload["constructed_theoretical_objects"]
        self.assertTrue(all(row["finite_partition_weight_witness"] for row in objs["gibbs_weight_rows"]))
        self.assertTrue(all(not row["physical_temperature_attached"] for row in objs["gibbs_weight_rows"]))
        self.assertTrue(all(row["detailed_balance_residual"] == 0.0 for row in objs["detailed_balance_transition_rows"]))
        self.assertTrue(all(not row["physical_bath_attached"] for row in objs["detailed_balance_transition_rows"]))
        self.assertTrue(all(row["kms_periodicity_proxy_witness"] for row in objs["kms_periodicity_proxy_rows"]))
        self.assertTrue(all(not row["physical_imaginary_time_clock_attached"] for row in objs["kms_periodicity_proxy_rows"]))
        self.assertTrue(all(row["finite_fdt_like_witness"] for row in objs["fluctuation_dissipation_proxy_rows"]))

    def test_negative_flags_and_recommendation(self):
        decision = self.payload["decision"]
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))
        self.assertIn("irreversibility/entropy-production", decision["next_honest_step"])
        self.assertIn("P3098/S2048", MD.read_text(encoding="utf-8"))
        self.assertIn("P3098/S2048", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3098/S2048", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3098", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
