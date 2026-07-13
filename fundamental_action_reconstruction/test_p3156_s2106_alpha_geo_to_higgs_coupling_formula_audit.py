import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3156_s2106_alpha_geo_to_higgs_coupling_formula_audit.py"
OUT = ROOT / "generated" / "p3156_s2106_alpha_geo_to_higgs_coupling_formula_audit.json"
MD = ROOT / "generated" / "p3156_s2106_alpha_geo_to_higgs_coupling_formula_audit.md"


class P3156AlphaGeoToHiggsCouplingFormulaAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_formula_counts_and_status(self):
        self.assertEqual(self.payload["status"], "P3156_ALPHA_GEO_TO_HIGGS_COUPLING_FORMULA_AUDIT_NO_STRICT_UNIT_SOURCE")
        counts = self.payload["finite_theorem"]["finite_counts"]
        self.assertEqual(counts["formula_schema_count"], 1)
        self.assertEqual(counts["formula_variant_rows"], 3)
        self.assertEqual(counts["dimensionally_valid_rows"], 2)
        self.assertEqual(counts["rows_importing_mass_scale"], 2)
        self.assertEqual(counts["accepted_strict_formula_rows"], 0)

    def test_dimension_obstruction(self):
        rows = {row["row"]: row for row in self.payload["formula_rows"]}
        self.assertFalse(rows["raw_dimensionless_alpha"]["dimensionally_valid_higgs_couplings"])
        self.assertTrue(rows["alpha_ratio_with_imported_mass_scale"]["dimensionally_valid_higgs_couplings"])
        self.assertTrue(rows["alpha_ratio_with_imported_mass_scale"]["imports_mass_scale"])
        self.assertFalse(any(row["accepted_strict_formula"] for row in self.payload["formula_rows"]))
        self.assertIn("Omega_dim/K_dim", self.payload["decision"]["next_honest_step"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_docs_updated(self):
        self.assertIn("P3156/S2106", MD.read_text(encoding="utf-8"))
        self.assertIn("P3156/S2106", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3156/S2106", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3156", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
