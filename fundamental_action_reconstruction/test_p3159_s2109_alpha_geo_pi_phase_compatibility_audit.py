import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3159_s2109_alpha_geo_pi_phase_compatibility_audit.py"
OUT = ROOT / "generated" / "p3159_s2109_alpha_geo_pi_phase_compatibility_audit.json"
MD = ROOT / "generated" / "p3159_s2109_alpha_geo_pi_phase_compatibility_audit.md"


class P3159AlphaGeoPiPhaseCompatibilityAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_exact_phase_section_and_status(self):
        self.assertEqual(self.payload["status"], "P3159_ALPHA_GEO_PI_PHASE_COMPATIBILITY_SCOPED_SECTION_NO_UNIT_SOURCE")
        cert = self.payload["numeric_certificate"]
        self.assertAlmostEqual(cert["alpha_geo"] * cert["A_phi"], cert["two_pi"], places=12)
        self.assertTrue(self.payload["decision"]["positive_scoped_flags"]["alpha_geo_pi_phase_compatibility_verified"])

    def test_no_closure_exports_and_rows(self):
        self.assertEqual(len(self.payload["phase_rows_n_1_to_12"]), 12)
        self.assertEqual(len(self.payload["continued_fraction_rows_alpha_over_2pi"]), 8)
        self.assertEqual(len(self.payload["compatibility_matrix"]), 6)
        self.assertEqual(self.payload["finite_theorem"]["accepted_strict_closures"], 0)
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))
        self.assertIn("Omega_M/K_dim", self.payload["decision"]["next_honest_step"])

    def test_docs_updated(self):
        self.assertIn("P3159/S2109", MD.read_text(encoding="utf-8"))
        self.assertIn("P3159/S2109", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3159/S2109", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3159", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
