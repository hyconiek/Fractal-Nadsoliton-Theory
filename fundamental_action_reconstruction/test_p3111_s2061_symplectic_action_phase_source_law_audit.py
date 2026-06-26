import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3111_s2061_symplectic_action_phase_source_law_audit.py"
OUT = ROOT / "generated" / "p3111_s2061_symplectic_action_phase_source_law_audit.json"
MD = ROOT / "generated" / "p3111_s2061_symplectic_action_phase_source_law_audit.md"


class P3111SymplecticActionPhaseSourceLawAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_predecessor_hash(self):
        self.assertEqual(self.payload["status"], "P3111_INTERNAL_PHASE_AREA_SECTION_POSITIVE_PHYSICAL_UNIT_NO_GO")
        self.assertTrue(self.payload["input_hashes"]["P3110"])

    def test_certificate_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 4)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["p3110_accepted_dimension_comparison_standards"], 0)
        self.assertEqual(cert["phase_periodicity_rows"], 8)
        self.assertEqual(cert["minimal_positive_phase_sections"], 1)
        self.assertEqual(cert["scale_orbit_rows"], 12)
        self.assertEqual(cert["candidate_source_laws"], 3)
        self.assertEqual(cert["candidate_gate_rows"], 24)
        self.assertEqual(cert["accepted_physical_unit_source_laws"], 0)

    def test_internal_phase_section_but_no_physical_unit(self):
        objs = self.payload["constructed_theoretical_objects"]
        phase_rows = objs["phase_periodicity_rows"]
        self.assertTrue(any(row["winding_n"] == 1 and row["internal_phase_section_selected"] for row in phase_rows))
        self.assertTrue(all(row["phase_residual_abs"] < 1e-12 for row in phase_rows))
        self.assertTrue(any(row["candidate"] == "minimal_internal_phase_periodicity_law" and row["fixes_unique_positive_phase_area"] for row in objs["candidate_source_law_rows"]))
        self.assertTrue(all(not row["accepted_physical_unit_source_law"] for row in objs["candidate_aggregate_certificate"]))

    def test_negative_flags_docs_and_recommendation(self):
        decision = self.payload["decision"]
        self.assertTrue(decision["positive_scoped_flags"]["minimal_internal_phase_area_section_exported"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))
        self.assertIn("internal calibration functional C_phi", decision["next_honest_step"])
        self.assertIn("P3111/S2061", MD.read_text(encoding="utf-8"))
        self.assertIn("P3111/S2061", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3111/S2061", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3111", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
