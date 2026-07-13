import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3152_s2102_y_sm_charge_unit_normalization_obstruction.py"
OUT = ROOT / "generated" / "p3152_s2102_y_sm_charge_unit_normalization_obstruction.json"
MD = ROOT / "generated" / "p3152_s2102_y_sm_charge_unit_normalization_obstruction.md"


class P3152YSmChargeUnitNormalizationObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_scale_orbit_is_invariant_but_not_unique(self):
        self.assertEqual(self.payload["status"], "P3152_YSM_CHARGE_UNIT_NORMALIZATION_SCALE_TORSOR_NO_STRICT_SOURCE")
        counts = self.payload["finite_theorem"]["finite_counts"]
        self.assertEqual(counts["sampled_nonzero_rational_scales"], 8)
        self.assertEqual(counts["invariant_scales"], 8)
        self.assertEqual(counts["distinct_Y_H_values"], 8)
        self.assertEqual(counts["strict_accepted_source_rows"], 0)

    def test_every_row_has_zero_residuals_and_varied_units(self):
        rows = self.payload["scale_orbit_rows"]
        self.assertTrue(all(row["all_yukawa_and_anomaly_residuals_zero"] for row in rows))
        self.assertIn("1/2", {row["Y_H"] for row in rows})
        self.assertIn("1", {row["Y_H"] for row in rows})
        self.assertGreater(len({row["electric_charge_gcd"] for row in rows}), 1)

    def test_no_strict_charge_source_and_docs_updated(self):
        self.assertFalse(any(row["fixes_scale"] and row["strict_source_law"] and row["noncircular"] for row in self.payload["candidate_source_rows"]))
        self.assertIn("P3153", self.payload["decision"]["next_honest_step"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))
        self.assertIn("P3152/S2102", MD.read_text(encoding="utf-8"))
        self.assertIn("P3152/S2102", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3152/S2102", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3152", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
