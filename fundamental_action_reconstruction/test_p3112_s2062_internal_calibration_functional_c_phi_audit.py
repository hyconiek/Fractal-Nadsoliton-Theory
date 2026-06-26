import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3112_s2062_internal_calibration_functional_c_phi_audit.py"
OUT = ROOT / "generated" / "p3112_s2062_internal_calibration_functional_c_phi_audit.json"
MD = ROOT / "generated" / "p3112_s2062_internal_calibration_functional_c_phi_audit.md"


class P3112InternalCalibrationFunctionalCPhiAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_predecessor_hash(self):
        self.assertEqual(self.payload["status"], "P3112_INTERNAL_C_PHI_CALIBRATION_FUNCTIONAL_BOUNDED_NO_GO")
        self.assertTrue(self.payload["input_hashes"]["P3111"])

    def test_certificate_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 4)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertTrue(cert["p3111_minimal_internal_phase_area_section_exported"])
        self.assertEqual(cert["candidate_C_phi_functionals"], 5)
        self.assertEqual(cert["scale_covariance_witness_rows"], 25)
        self.assertEqual(cert["action_length_time_induction_rows"], 15)
        self.assertEqual(cert["candidate_gate_rows"], 45)
        self.assertEqual(cert["accepted_internal_dimensionful_calibration_functionals"], 0)

    def test_c_phi_candidates_and_no_go(self):
        objs = self.payload["constructed_theoretical_objects"]
        candidates = objs["candidate_C_phi_functionals"]
        self.assertTrue(any(row["candidate"] == "alpha_geo_phase_normalized_C_phi" and row["value_at_A_phi"] == 1.0 for row in candidates))
        self.assertTrue(any(row["candidate"] == "imported_hbar_planck_C_phi" and not row["standard_physics_import_free"] for row in candidates))
        self.assertTrue(all(not row["accepted_internal_dimensionful_calibration_functional"] for row in objs["candidate_aggregate_certificate"]))
        self.assertTrue(all(not row["target_calibration_induced"] for row in objs["action_length_time_induction_rows"]))

    def test_negative_flags_docs_and_recommendation(self):
        decision = self.payload["decision"]
        self.assertTrue(decision["positive_scoped_flags"]["candidate_C_phi_functionals_constructed"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))
        self.assertIn("dimensionful reference-carrier source law U_action", decision["next_honest_step"])
        self.assertIn("P3112/S2062", MD.read_text(encoding="utf-8"))
        self.assertIn("P3112/S2062", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3112/S2062", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3112", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
