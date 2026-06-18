import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
P2885_SCRIPT = ROOT / "p2885_s1835_invariant_quadratic_action_9_over_5_selection_no_go_audit.py"
SCRIPT = ROOT / "p2886_s1836_external_unit_measure_action_density_inventory_no_go_audit.py"
JSON_PATH = ROOT / "generated" / "p2886_s1836_external_unit_measure_action_density_inventory_no_go_audit.json"
MD_PATH = ROOT / "generated" / "p2886_s1836_external_unit_measure_action_density_inventory_no_go_audit.md"


class P2886ExternalUnitMeasureActionDensityInventoryNoGoAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for script in (P2885_SCRIPT, SCRIPT, SCRIPT):
            subprocess.run([sys.executable, str(script)], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["external_unit_measure_action_density_inventory_no_go_audit"]
        cls.facts = cls.payload["acceptance_matrix"]["facts"]
        cls.flags = cls.payload["decision"]["negative_export_flags"]

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2886_EXTERNAL_UNIT_MEASURE_ACTION_DENSITY_INVENTORY_NO_GO_AUDIT_NO_CLOSURE")
        self.assertEqual(self.audit["input_status_rechecked"], "P2885_INVARIANT_QUADRATIC_ACTION_9_OVER_5_SELECTION_NO_GO_AUDIT_NO_CLOSURE")
        self.assertTrue(self.facts["p2885_rechecked"])

    def test_inventory_is_nonempty_and_relevant(self):
        self.assertGreater(self.audit["generated_json_file_count"], 0)
        self.assertGreater(self.audit["relevant_record_count"], 0)
        self.assertTrue(self.facts["generated_json_inventory_nonempty"])
        self.assertTrue(self.facts["unit_measure_or_action_density_terms_found"])

    def test_no_positive_external_export_found(self):
        self.assertEqual(self.audit["positive_export_record_count"], 0)
        self.assertTrue(self.facts["no_positive_external_unit_measure_or_action_density_export_found"])
        self.assertTrue(self.facts["p2885_obligation_remains_unsatisfied"])

    def test_false_export_flags(self):
        self.assertFalse(any(self.flags.values()))
        self.assertFalse(self.payload["acceptance_matrix"]["accepted_as_external_unit_measure_or_action_density_export"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_unit_bearing_action_density"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_variational_chain_rule_to_ltotal"])

    def test_documents_updated(self):
        self.assertIn("P2886/S1836", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2886/S1836", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2886/S1836", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2886", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
