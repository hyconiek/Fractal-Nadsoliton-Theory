import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2957_s1907_positive_beta_scale_unit_source_obstruction.py"
OUT = ROOT / "generated" / "p2957_s1907_positive_beta_scale_unit_source_obstruction.json"
MD = ROOT / "generated" / "p2957_s1907_positive_beta_scale_unit_source_obstruction.md"


class P2957PositiveBetaScaleUnitSourceObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input_hashes(self):
        self.assertEqual(self.payload["status"], "P2957_POSITIVE_BETA_SCALE_UNIT_SOURCE_OBSTRUCTION_NO_STRICT_EXPORT")
        for key in ["P2692", "P2950", "P2951", "P2956"]:
            self.assertIsNotNone(self.payload["input_hashes"][key])

    def test_scale_covariance_rows(self):
        cert = self.payload["positive_beta_scale_unit_certificate"]
        self.assertEqual(cert["eta"]["as_string"], "9/5")
        self.assertEqual(cert["sample_positive_beta_representative_count"], 5)
        self.assertTrue(cert["all_rows_positive_and_covariant"])
        self.assertFalse(cert["any_row_selects_canonical_unit"])
        rows = self.payload["constructed_theoretical_objects"]["scale_covariance_rows"]
        self.assertTrue(all(row["orbit_covariance_available"] for row in rows))
        self.assertFalse(any(row["canonical_unit_selected_by_row"] for row in rows))

    def test_obligations_and_acceptance_matrix(self):
        obligations = {row["obligation"]: row for row in self.payload["constructed_theoretical_objects"]["source_obligation_rows"]}
        self.assertTrue(obligations["exact_ratio_package_eta_available"]["satisfied"])
        self.assertTrue(obligations["positive_beta_orbit_covariance_available"]["satisfied"])
        self.assertTrue(obligations["p2951_positive_beta_scale_unit_atom_named"]["satisfied"])
        self.assertFalse(obligations["target_independent_positive_beta_source_exported"]["satisfied"])
        self.assertFalse(obligations["canonical_length_uv_unit_exported"]["satisfied"])
        self.assertFalse(obligations["unit_bearing_nonproxy_coupling_exported"]["satisfied"])
        self.assertFalse(obligations["scale_orbit_quotient_selects_unique_positive_unit"]["satisfied"])
        matrix = self.payload["constructed_theoretical_objects"]["finite_acceptance_matrix"]
        self.assertEqual(len(matrix), 64)
        self.assertEqual(sum(1 for row in matrix if row["accepts_strict_positive_beta_scale_unit_source"]), 1)

    def test_certificate_nonpromotion_and_docs(self):
        cert = self.payload["positive_beta_scale_unit_certificate"]
        self.assertFalse(cert["target_independent_positive_beta_source_exported"])
        self.assertFalse(cert["canonical_length_uv_unit_exported"])
        self.assertFalse(cert["unit_bearing_nonproxy_coupling_exported"])
        self.assertFalse(cert["scale_orbit_quotient_selects_unique_positive_unit"])
        self.assertFalse(cert["p2951_positive_beta_scale_unit_atom_discharged"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertIn("P2957/S1907", MD.read_text(encoding="utf-8"))
        self.assertIn("P2957/S1907", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2957/S1907", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2957", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
