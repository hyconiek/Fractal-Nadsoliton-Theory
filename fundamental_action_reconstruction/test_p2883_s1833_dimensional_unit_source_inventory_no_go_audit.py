import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
P2882_SCRIPT = ROOT / "p2882_s1832_euler_source_ratio_9_over_5_forcing_no_go_audit.py"
SCRIPT = ROOT / "p2883_s1833_dimensional_unit_source_inventory_no_go_audit.py"
JSON_PATH = ROOT / "generated" / "p2883_s1833_dimensional_unit_source_inventory_no_go_audit.json"
MD_PATH = ROOT / "generated" / "p2883_s1833_dimensional_unit_source_inventory_no_go_audit.md"


class P2883DimensionalUnitSourceInventoryNoGoAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for script in (P2882_SCRIPT, SCRIPT):
            subprocess.run([sys.executable, str(script)], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["dimensional_unit_source_inventory_no_go_audit"]
        cls.facts = cls.payload["acceptance_matrix"]["facts"]
        cls.flags = cls.payload["decision"]["negative_export_flags"]

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2883_DIMENSIONAL_UNIT_SOURCE_INVENTORY_NO_GO_AUDIT_NO_CLOSURE")
        self.assertEqual(self.audit["input_status_rechecked"], "P2882_EULER_SOURCE_RATIO_9_OVER_5_FORCING_NO_GO_AUDIT_NO_CLOSURE")

    def test_inventory_is_nonempty_and_relevant(self):
        self.assertGreater(self.audit["generated_json_file_count"], 0)
        self.assertGreater(self.audit["relevant_record_count"], 0)
        self.assertTrue(self.facts["generated_json_inventory_nonempty"])
        self.assertTrue(self.facts["relevant_inventory_records_found"])

    def test_positive_hits_are_insufficient(self):
        self.assertGreater(self.audit["positive_relevant_record_count"], 0)
        self.assertIn("eta/damping/unit-premise", self.audit["positive_hits_classification"])
        self.assertEqual(self.audit["positive_source_stiffness_9_to_5_record_count"], 0)
        self.assertTrue(self.facts["no_positive_source_stiffness_9_to_5_export_found"])
        self.assertTrue(self.facts["positive_9_over_5_eta_hits_do_not_export_source_stiffness_law"])

    def test_false_export_flags(self):
        self.assertFalse(any(self.flags.values()))
        self.assertFalse(self.payload["acceptance_matrix"]["accepted_as_current_dimensional_unit_source_export"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_source_stiffness_ratio_9_to_5"])

    def test_documents_updated(self):
        self.assertIn("P2883/S1833", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2883/S1833", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2883/S1833", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2883", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
