import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2750_s1700_concrete_odd_source_sign_value_inventory_no_go.py"
OUT = ROOT / "generated" / "p2750_s1700_concrete_odd_source_sign_value_inventory_no_go.json"
MD = ROOT / "generated" / "p2750_s1700_concrete_odd_source_sign_value_inventory_no_go.md"


class P2750ConcreteOddSourceInventoryNoGoTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.scan = cls.payload["content_evidence_scan"]
        cls.inventory = cls.payload["candidate_inventory"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_content_first_grep_boundaries_present(self):
        self.assertEqual(self.payload["status"], "P2750_CONCRETE_ODD_SOURCE_SIGN_VALUE_INVENTORY_NO_GO")
        self.assertTrue(self.scan["all_patterns_have_hits"])
        for lane in ("post_p2749_missing_premise", "candidate_signed_observables", "p2721_boundary", "closure_forbidden"):
            self.assertGreater(self.scan["hit_counts"][lane], 0)

    def test_inventory_loads_signed_candidates_but_accepts_none(self):
        self.assertEqual(self.inventory["candidate_file_count"], 10)
        self.assertGreaterEqual(self.inventory["existing_candidate_file_count"], 9)
        self.assertEqual(self.inventory["accepted_candidate_count"], 0)
        self.assertEqual(self.inventory["accepted_artifacts"], [])
        self.assertTrue(any(row["has_candidate_signed_object"] for row in self.inventory["rows"]))
        self.assertTrue(any(row["has_concrete_numeric_signal"] for row in self.inventory["rows"]))

    def test_acceptance_matrix_preserves_p2749_blocker(self):
        self.assertFalse(self.acceptance["accepted_as_p2749_completion"])
        self.assertFalse(self.acceptance["facts"]["accepted_concrete_strict_source_sign_value_found"])
        self.assertFalse(self.acceptance["facts"]["unique_p2721_coupling_polarity_theorem_found"])
        self.assertIn("accepted_concrete_strict_source_sign_value_found", self.acceptance["missing_criteria"])
        self.assertIn("unique_p2721_coupling_polarity_theorem_found", self.acceptance["missing_criteria"])

    def test_negative_export_flags_remain_false(self):
        flags = self.payload["decision"]["negative_export_flags"]
        self.assertTrue(flags)
        self.assertTrue(all(value is False for value in flags.values()))

    def test_recommendation_and_docs_are_updated(self):
        self.assertIn("P2697-P2750", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2750/S1700", MD.read_text(encoding="utf-8"))
        self.assertIn("P2750/S1700", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2750/S1700", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2750", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
