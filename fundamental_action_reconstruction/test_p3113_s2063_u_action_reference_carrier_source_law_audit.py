import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3113_s2063_u_action_reference_carrier_source_law_audit.py"
OUT = ROOT / "generated" / "p3113_s2063_u_action_reference_carrier_source_law_audit.json"
MD = ROOT / "generated" / "p3113_s2063_u_action_reference_carrier_source_law_audit.md"


class P3113UActionReferenceCarrierSourceLawAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_predecessor_hash(self):
        self.assertEqual(self.payload["status"], "P3113_U_ACTION_REFERENCE_CARRIER_SOURCE_LAW_BOUNDED_NO_GO")
        self.assertTrue(self.payload["input_hashes"]["P3112"])

    def test_certificate_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 4)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["p3112_accepted_internal_dimensionful_calibration_functionals"], 0)
        self.assertEqual(cert["candidate_U_action_source_laws"], 6)
        self.assertEqual(cert["scale_orbit_section_rows"], 30)
        self.assertEqual(cert["dimensional_balance_rows"], 18)
        self.assertEqual(cert["C_phi_coupling_rows"], 6)
        self.assertEqual(cert["candidate_gate_rows"], 54)
        self.assertEqual(cert["accepted_U_action_source_laws"], 0)

    def test_candidate_matrix_and_rejections(self):
        objs = self.payload["constructed_theoretical_objects"]
        candidates = objs["candidate_U_action_source_laws"]
        self.assertTrue(any(row["candidate"] == "formal_declared_U_action_symbol" and row["nonzero_action_dimension_exported"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "imported_planck_action_carrier" and not row["standard_physics_import_free"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "pure_phase_reference_carrier" and row["C_phi_A_phi_coupling_theorem_exported"] for row in candidates))
        self.assertTrue(all(not row["accepted_U_action_source_law"] for row in objs["candidate_aggregate_certificate"]))
        self.assertTrue(all(not row["coupling_accepted"] for row in objs["C_phi_coupling_rows"]))

    def test_negative_flags_docs_and_recommendation(self):
        decision = self.payload["decision"]
        self.assertTrue(decision["positive_scoped_flags"]["candidate_U_action_source_laws_constructed"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))
        self.assertIn("dimensional triad source law D_phi", decision["next_honest_step"])
        self.assertIn("P3113/S2063", MD.read_text(encoding="utf-8"))
        self.assertIn("P3113/S2063", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3113/S2063", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3113", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
