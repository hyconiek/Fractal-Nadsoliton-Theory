import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
P2884_SCRIPT = ROOT / "p2884_s1834_finite_invariant_ratio_9_to_5_source_law_no_go_audit.py"
SCRIPT = ROOT / "p2885_s1835_invariant_quadratic_action_9_over_5_selection_no_go_audit.py"
JSON_PATH = ROOT / "generated" / "p2885_s1835_invariant_quadratic_action_9_over_5_selection_no_go_audit.json"
MD_PATH = ROOT / "generated" / "p2885_s1835_invariant_quadratic_action_9_over_5_selection_no_go_audit.md"


class P2885InvariantQuadraticActionSelectionNoGoAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for script in (P2884_SCRIPT, SCRIPT):
            subprocess.run([sys.executable, str(script)], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["invariant_quadratic_action_9_over_5_selection_no_go_audit"]
        cls.facts = cls.payload["acceptance_matrix"]["facts"]
        cls.flags = cls.payload["decision"]["negative_export_flags"]

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2885_INVARIANT_QUADRATIC_ACTION_9_OVER_5_SELECTION_NO_GO_AUDIT_NO_CLOSURE")
        self.assertEqual(self.audit["input_status_rechecked"], "P2884_FINITE_INVARIANT_RATIO_9_TO_5_SOURCE_LAW_NO_GO_AUDIT_NO_CLOSURE")
        self.assertTrue(self.facts["p2884_rechecked"])

    def test_quadratic_action_family_is_complete_over_p2884_expressions(self):
        self.assertGreater(self.audit["candidate_action_count"], 0)
        self.assertEqual(self.audit["candidate_action_count"], 446224)
        self.assertTrue(self.facts["quadratic_action_family_nonempty"])

    def test_target_actions_import_ratio_and_are_not_unique(self):
        self.assertGreater(self.audit["target_action_record_count"], 1)
        self.assertEqual(self.audit["target_without_primitive_9_to_5_ratio_count"], 0)
        self.assertGreater(self.audit["distinct_target_expression_pair_count"], 1)
        self.assertTrue(self.facts["every_target_imports_primitive_9_to_5_ratio"])
        self.assertTrue(self.facts["target_action_not_unique"])
        self.assertTrue(self.facts["no_exported_unit_weight_or_measure_selects_target_action"])

    def test_false_export_flags(self):
        self.assertFalse(any(self.flags.values()))
        self.assertFalse(self.payload["acceptance_matrix"]["accepted_as_unit_dimensional_action_functional_selecting_9_over_5"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_unique_source_stiffness_ratio_9_to_5"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_variational_chain_rule_to_ltotal"])

    def test_documents_updated(self):
        self.assertIn("P2885/S1835", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2885/S1835", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2885/S1835", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2885", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
