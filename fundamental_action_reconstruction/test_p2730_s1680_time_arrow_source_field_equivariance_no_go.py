import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
JSON_PATH = GEN / "p2730_s1680_time_arrow_source_field_equivariance_no_go.json"
MD_PATH = GEN / "p2730_s1680_time_arrow_source_field_equivariance_no_go.md"
SCRIPT = ROOT / "p2730_s1680_time_arrow_source_field_equivariance_no_go.py"


class P2730TimeArrowSourceFieldEquivarianceNoGoTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["time_arrow_source_field_audit"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_exhausts_all_z12_tau_fields(self):
        self.assertEqual(self.payload["status"], "P2730_TIME_ARROW_SOURCE_FIELD_EQUIVARIANCE_NO_GO")
        self.assertEqual(self.audit["field_count"], 4096)
        self.assertEqual(self.audit["source_orbits_under_aut_z12"], [[0], [1, 5, 7, 11], [2, 10], [3, 9], [4, 8], [6]])

    def test_translation_safe_fields_are_only_paired_constants(self):
        self.assertEqual(self.audit["translation_invariant_field_count"], 2)
        self.assertEqual(self.audit["translation_and_aut_invariant_field_count"], 2)
        self.assertEqual(self.audit["translation_and_aut_invariant_fields"], [[-1] * 12, [1] * 12])
        self.assertEqual(self.audit["sign_unpaired_translation_invariant_field_count"], 0)

    def test_aut_invariant_fields_are_sign_paired(self):
        self.assertEqual(self.audit["aut_invariant_field_count"], 64)
        self.assertEqual(self.audit["nonconstant_aut_invariant_field_count"], 62)
        self.assertEqual(self.audit["sign_unpaired_aut_invariant_field_count"], 0)

    def test_acceptance_blocks_strict_source_and_closure_flags(self):
        self.assertFalse(self.acceptance["accepted_as_strict_time_arrow_source"])
        self.assertIn("translation_gauge_safe_tau_sign_unpaired", self.acceptance["missing_criteria"])
        self.assertIn("strict_nonpremise_tau_sign_selected", self.acceptance["missing_criteria"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_documentation_updated(self):
        self.assertIn("P2730/S1680", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2730/S1680", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2730/S1680", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2730", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
