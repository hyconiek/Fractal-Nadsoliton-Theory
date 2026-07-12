import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3133_s2083_legacy_beta_tors_helical_lock_defect_audit.py"
OUT = ROOT / "generated" / "p3133_s2083_legacy_beta_tors_helical_lock_defect_audit.json"
MD = ROOT / "generated" / "p3133_s2083_legacy_beta_tors_helical_lock_defect_audit.md"


class P3133LegacyBetaTorsHelicalLockDefectAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P3133_LEGACY_BETA_TORS_HELICAL_LOCK_DEFECT_GAP_CERTIFICATE")
        self.assertTrue(self.payload["input_hashes"]["P3132"])
        self.assertTrue(self.payload["input_hashes"]["K1"])
        self.assertTrue(self.payload["input_hashes"]["QW1729"])

    def test_legacy_rows_capture_gap_without_export(self):
        rows = self.payload["comparison_rows"]
        self.assertEqual(len(rows), 4)
        sinusoid = next(row for row in rows if row["legacy_atom"] == "cos(omega*d+phi) sinusoid")
        beta = next(row for row in rows if row["legacy_atom"] == "beta_tors linear denominator")
        self.assertTrue(sinusoid["D_HL_nontranslation_defect"])
        self.assertFalse(sinusoid["nonzero_inversion_odd_value"])
        self.assertFalse(beta["D_HL_nontranslation_defect"])
        self.assertIn("scalar damping/torsion marker", beta["fresh_gap"])

    def test_certificate_and_recommendation(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["legacy_atoms_tested"], 4)
        self.assertEqual(cert["accepted_D_HL_sources"], 0)
        self.assertTrue(cert["sinusoid_retains_phase_trace"])
        self.assertTrue(cert["beta_tors_retains_scalar_torsion_damping_trace"])
        self.assertFalse(cert["legacy_role_transfer_started"])
        decision = self.payload["decision"]
        self.assertTrue(decision["positive_scoped_flags"]["strict_operational_absorption_gap_identified"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))
        self.assertIn("explicit candidate formula for D_HL", decision["next_honest_step"])

    def test_docs_updated(self):
        self.assertIn("P3133/S2083", MD.read_text(encoding="utf-8"))
        self.assertIn("P3133/S2083", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3133/S2083", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3133", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
