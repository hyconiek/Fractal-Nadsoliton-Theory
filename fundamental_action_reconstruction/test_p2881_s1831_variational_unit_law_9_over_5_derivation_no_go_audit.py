import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
P2880_SCRIPT = ROOT / "p2880_s1830_origin_pinned_9_over_5_coefficient_import_no_go_audit.py"
SCRIPT = ROOT / "p2881_s1831_variational_unit_law_9_over_5_derivation_no_go_audit.py"
JSON_PATH = ROOT / "generated" / "p2881_s1831_variational_unit_law_9_over_5_derivation_no_go_audit.json"
MD_PATH = ROOT / "generated" / "p2881_s1831_variational_unit_law_9_over_5_derivation_no_go_audit.md"


class P2881VariationalUnitLawNineOverFiveDerivationNoGoAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for script in (P2880_SCRIPT, SCRIPT):
            subprocess.run([sys.executable, str(script)], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["variational_unit_law_9_over_5_derivation_no_go_audit"]
        cls.facts = cls.payload["acceptance_matrix"]["facts"]
        cls.flags = cls.payload["decision"]["negative_export_flags"]

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2881_VARIATIONAL_UNIT_LAW_9_OVER_5_DERIVATION_NO_GO_AUDIT_NO_CLOSURE")
        self.assertEqual(self.audit["input_status_rechecked"], "P2880_ORIGIN_PINNED_9_OVER_5_COEFFICIENT_IMPORT_NO_GO_AUDIT_NO_CLOSURE")

    def test_candidate_grid_and_family_size(self):
        self.assertEqual(self.audit["stencil_count"], 19**3)
        self.assertEqual(self.audit["constraint_count"], 15)
        self.assertEqual(self.audit["objective_count"], 5)
        self.assertEqual(self.audit["candidate_count"], 75)
        self.assertTrue(self.facts["candidate_count_checked"])

    def test_no_unique_variational_derivation_of_9_over_5(self):
        self.assertEqual(self.audit["derives_center_9_over_5_candidate_count"], 0)
        self.assertTrue(self.facts["no_candidate_derives_unique_center_9_over_5"])
        self.assertGreater(self.audit["unique_minimizer_candidate_count"], 0)
        self.assertTrue(self.facts["some_unique_minimizers_exist_but_not_9_over_5"])

    def test_center_9_over_5_is_not_accepted_when_seen(self):
        self.assertGreaterEqual(self.audit["center_9_over_5_nonunique_minimum_candidate_count"], 0)
        self.assertTrue(self.facts["center_9_over_5_only_appears_in_nonunique_minima_when_seen"])
        self.assertFalse(self.payload["acceptance_matrix"]["accepted_as_independent_variational_unit_law_for_9_over_5"])

    def test_false_export_flags(self):
        self.assertFalse(any(self.flags.values()))
        self.assertFalse(self.payload["acceptance_matrix"]["exports_unit_bearing_9_over_5_coupling_theorem"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_unit_bearing_action_density"])

    def test_documents_updated(self):
        self.assertIn("P2881/S1831", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2881/S1831", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2881/S1831", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2881", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
