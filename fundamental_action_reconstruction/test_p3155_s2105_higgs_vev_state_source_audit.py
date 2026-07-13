import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3155_s2105_higgs_vev_state_source_audit.py"
OUT = ROOT / "generated" / "p3155_s2105_higgs_vev_state_source_audit.json"
MD = ROOT / "generated" / "p3155_s2105_higgs_vev_state_source_audit.md"


class P3155HiggsVevStateSourceAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_stationary_structure_and_counts(self):
        self.assertEqual(self.payload["status"], "P3155_HIGGS_VEV_STATE_SOURCE_AUDIT_NO_STRICT_SOURCE")
        counts = self.payload["finite_theorem"]["finite_counts"]
        self.assertEqual(counts["stationary_points"], 3)
        self.assertEqual(counts["nonzero_stationary_points"], 2)
        self.assertEqual(counts["scalar_candidate_rows"], 5)
        self.assertEqual(counts["positive_dimensionless_candidate_rows"], 3)
        self.assertEqual(counts["accepted_strict_vev_source_rows"], 0)
        self.assertEqual(self.payload["quartic_potential"]["nonzero_vev_formula"], "v^2 = mu2/lambda")

    def test_no_candidate_exports_full_vev_source(self):
        candidates = self.payload["scalar_candidate_rows"]
        self.assertFalse(any(row["fixes_mu2"] and row["fixes_lambda"] and row["fixes_ratio_v2"] and row["fixes_dimensionful_vev_unit"] and row["noncircular_state_equation"] for row in candidates))
        self.assertIn("below physical model status", self.payload["decision"]["model_assessment"])
        self.assertIn("P3156", self.payload["decision"]["next_honest_step"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_docs_updated(self):
        self.assertIn("P3155/S2105", MD.read_text(encoding="utf-8"))
        self.assertIn("P3155/S2105", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3155/S2105", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3155", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
