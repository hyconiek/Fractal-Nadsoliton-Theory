import json, subprocess, sys, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3100_s2050_open_system_bath_preparation_source_obstruction_audit.py"
OUT = ROOT / "generated" / "p3100_s2050_open_system_bath_preparation_source_obstruction_audit.json"
MD = ROOT / "generated" / "p3100_s2050_open_system_bath_preparation_source_obstruction_audit.md"

class P3100OpenSystemBathPreparationSourceObstructionAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_predecessor_hash(self):
        self.assertEqual(self.payload["status"], "P3100_OPEN_SYSTEM_BATH_PREPARATION_SOURCE_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertTrue(self.payload["input_hashes"]["P3099"])

    def test_certificate_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["p3099_accepted_nonimported_thermalization_sources"], 0)
        self.assertEqual(cert["relative_entropy_monotonicity_rows"], 96)
        self.assertEqual(cert["relative_entropy_rows_monotone"], 96)
        self.assertEqual(cert["entropy_production_proxy_rows"], 528)
        self.assertEqual(cert["stochastic_semigroup_proxy_rows"], 4)
        self.assertEqual(cert["modal_flux_relaxation_rows"], 24)
        self.assertEqual(cert["bath_preparation_candidates"], 5)
        self.assertEqual(cert["required_gates"], 8)
        self.assertEqual(cert["candidate_gate_rows"], 40)
        self.assertEqual(cert["accepted_nonimported_bath_preparation_sources"], 0)

    def test_witness_rows_are_formal_and_unsourced(self):
        objs = self.payload["constructed_theoretical_objects"]
        self.assertTrue(all(row["relative_entropy_monotonicity_witness"] for row in objs["relative_entropy_monotonicity_rows"]))
        self.assertTrue(all(row["monotone_nonincreasing_from_previous"] for row in objs["relative_entropy_monotonicity_rows"]))
        self.assertTrue(all(not row["physical_time_step_attached"] for row in objs["relative_entropy_monotonicity_rows"]))
        self.assertTrue(all(row["finite_entropy_production_proxy_witness"] for row in objs["entropy_production_proxy_rows"]))
        self.assertTrue(all(not row["physical_bath_attached"] for row in objs["entropy_production_proxy_rows"]))
        self.assertTrue(all(row["finite_stochastic_semigroup_proxy_witness"] for row in objs["stochastic_semigroup_proxy_rows"]))
        self.assertTrue(all(row["modal_relaxation_witness"] for row in objs["modal_flux_relaxation_rows"]))

    def test_negative_flags_and_recommendation(self):
        decision = self.payload["decision"]
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))
        self.assertIn("detector/readout calibration", decision["next_honest_step"])
        self.assertIn("P3100/S2050", MD.read_text(encoding="utf-8"))
        self.assertIn("P3100/S2050", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3100/S2050", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3100", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
