import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3157_s2107_omega_dim_mass_unit_torsor_audit.py"
OUT = ROOT / "generated" / "p3157_s2107_omega_dim_mass_unit_torsor_audit.json"
MD = ROOT / "generated" / "p3157_s2107_omega_dim_mass_unit_torsor_audit.md"


class P3157OmegaDimMassUnitTorsorAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_torsor_counts_and_status(self):
        self.assertEqual(self.payload["status"], "P3157_OMEGA_DIM_MASS_UNIT_TORSOR_FORMAL_CARRIER_NO_STRICT_SOURCE")
        counts = self.payload["finite_theorem"]["finite_counts"]
        self.assertEqual(counts["torsor_scale_rows"], 4)
        self.assertEqual(counts["dimensionally_valid_rows"], 4)
        self.assertEqual(counts["canonical_rows_by_gauge_label"], 1)
        self.assertEqual(counts["strict_source_selected_rows"], 0)
        self.assertEqual(counts["gates_satisfied"], 2)

    def test_scale_orbit_not_sourced(self):
        rows = self.payload["torsor_rows"]
        self.assertTrue(all(row["dimensionally_valid"] for row in rows))
        self.assertTrue(all(row["same_dimensionless_ratio"] for row in rows))
        self.assertFalse(any(row["strict_source_law_for_c"] for row in rows))
        self.assertIn("P3116", self.payload["decision"]["next_honest_step"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_docs_updated(self):
        self.assertIn("P3157/S2107", MD.read_text(encoding="utf-8"))
        self.assertIn("P3157/S2107", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3157/S2107", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3157", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
