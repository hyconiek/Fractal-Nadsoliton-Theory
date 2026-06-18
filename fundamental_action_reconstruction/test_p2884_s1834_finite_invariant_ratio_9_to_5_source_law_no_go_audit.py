import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
P2883_SCRIPT = ROOT / "p2883_s1833_dimensional_unit_source_inventory_no_go_audit.py"
SCRIPT = ROOT / "p2884_s1834_finite_invariant_ratio_9_to_5_source_law_no_go_audit.py"
JSON_PATH = ROOT / "generated" / "p2884_s1834_finite_invariant_ratio_9_to_5_source_law_no_go_audit.json"
MD_PATH = ROOT / "generated" / "p2884_s1834_finite_invariant_ratio_9_to_5_source_law_no_go_audit.md"


class P2884FiniteInvariantRatioSourceLawNoGoAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for script in (P2883_SCRIPT, SCRIPT):
            subprocess.run([sys.executable, str(script)], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["finite_invariant_ratio_9_to_5_source_law_no_go_audit"]
        cls.facts = cls.payload["acceptance_matrix"]["facts"]
        cls.flags = cls.payload["decision"]["negative_export_flags"]

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2884_FINITE_INVARIANT_RATIO_9_TO_5_SOURCE_LAW_NO_GO_AUDIT_NO_CLOSURE")
        self.assertEqual(self.audit["input_status_rechecked"], "P2883_DIMENSIONAL_UNIT_SOURCE_INVENTORY_NO_GO_AUDIT_NO_CLOSURE")
        self.assertTrue(self.facts["p2883_rechecked"])

    def test_finite_expression_and_ratio_family(self):
        self.assertGreater(self.audit["expression_count"], 0)
        self.assertEqual(self.audit["ratio_candidate_count"], self.audit["expression_count"] ** 2)
        self.assertTrue(self.facts["finite_invariant_expression_family_nonempty"])

    def test_9_to_5_is_representable_but_not_selected(self):
        self.assertGreater(self.audit["target_ratio_record_count"], 1)
        self.assertTrue(self.facts["target_ratio_representable_in_invariant_algebra"])
        self.assertTrue(self.facts["target_ratio_not_unique_or_selected"])
        self.assertTrue(self.facts["no_exported_selector_for_target_pair"])
        self.assertTrue(self.facts["no_unit_dimension_attached_to_ratio"])

    def test_false_export_flags(self):
        self.assertFalse(any(self.flags.values()))
        self.assertFalse(self.payload["acceptance_matrix"]["accepted_as_strict_invariant_source_law_for_9_to_5"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_source_stiffness_ratio_9_to_5"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_unit_bearing_action_density"])

    def test_documents_updated(self):
        self.assertIn("P2884/S1834", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2884/S1834", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2884/S1834", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2884", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
