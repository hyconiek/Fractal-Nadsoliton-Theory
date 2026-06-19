import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2900_s1850_defect_placement_source_law_inventory_no_go.py"
JSON_PATH = ROOT / "generated" / "p2900_s1850_defect_placement_source_law_inventory_no_go.json"
MD_PATH = ROOT / "generated" / "p2900_s1850_defect_placement_source_law_inventory_no_go.md"


class P2900DefectPlacementSourceLawInventoryNoGoTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.scan = cls.payload["defect_placement_source_law_inventory"]["scan"]
        cls.acceptance = cls.payload["acceptance_matrix"]
        cls.flags = cls.payload["decision"]["negative_export_flags"]

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2900_DEFECT_PLACEMENT_SOURCE_LAW_INVENTORY_NO_GO_NO_CLOSURE")
        self.assertTrue(self.acceptance["p2899_rechecked"])
        self.assertEqual(
            self.payload["defect_placement_source_law_inventory"]["input_status_rechecked"],
            "P2899_POST_DEFECT_POTENTIAL_READINESS_MATRIX_NO_CLOSURE",
        )

    def test_no_coupled_positive_export(self):
        self.assertGreater(self.scan["generated_json_file_count"], 0)
        self.assertGreater(self.scan["term_occurrence_count"], 0)
        self.assertEqual(self.scan["coupled_source_and_coupling_artifact_count"], 0)
        self.assertEqual(self.scan["coupled_source_coupling_and_closure_artifact_count"], 0)
        self.assertFalse(self.acceptance["exact_missing_object_found"])
        self.assertTrue(self.acceptance["accepted_as_current_artifact_inventory_no_go"])

    def test_false_export_flags(self):
        self.assertFalse(any(self.flags.values()))
        self.assertFalse(self.payload["decision"]["negative_export_flags"]["toe_closure_exported"])
        self.assertFalse(self.payload["decision"]["negative_export_flags"]["ltotal_exported"])

    def test_documents_updated(self):
        self.assertIn("P2900/S1850", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2900/S1850", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2900/S1850", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2900", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
