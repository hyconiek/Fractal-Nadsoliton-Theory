import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2748_s1698_absence_of_selector_self_synchronization_no_go.py"
OUT = ROOT / "generated" / "p2748_s1698_absence_of_selector_self_synchronization_no_go.json"
MD = ROOT / "generated" / "p2748_s1698_absence_of_selector_self_synchronization_no_go.md"


class P2748AbsenceOfSelectorSelfSynchronizationNoGoTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.scan = cls.payload["content_evidence_scan"]
        cls.audit = cls.payload["synchronization_audit"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_content_first_grep_detects_absence_selector_boundaries(self):
        self.assertEqual(self.payload["status"], "P2748_ABSENCE_OF_SELECTOR_SELF_SYNCHRONIZATION_NO_GO")
        self.assertEqual(self.scan["content_pattern_count"], 4)
        self.assertTrue(self.scan["all_patterns_have_hits"])
        self.assertGreater(self.scan["hit_counts"]["post_p2747_boundary"], 0)
        self.assertGreater(self.scan["hit_counts"]["absence_selector_language"], 0)

    def test_absence_data_have_no_equivariant_selector_map(self):
        self.assertEqual(self.audit["orientation_reversing_units"], [7, 11])
        self.assertEqual(self.audit["singleton_absence_map_count"], 2)
        self.assertEqual(self.audit["singleton_absence_equivariant_map_count"], 0)
        self.assertEqual(self.audit["bit_absence_map_count"], 4)
        self.assertEqual(self.audit["bit_absence_equivariant_map_count"], 0)

    def test_odd_signed_source_is_explicitly_extra_structure(self):
        self.assertEqual(self.audit["odd_signed_source_equivariant_map_count"], 2)
        self.assertEqual(len(self.audit["odd_signed_source_equivariant_maps"]), 2)
        self.assertIn("not information about absence alone", self.audit["finite_theorem"])

    def test_acceptance_blocks_self_synchronization_export(self):
        self.assertFalse(self.acceptance["accepted_as_self_synchronizing_selector"])
        self.assertIn("absence_singleton_exports_equivariant_selector", self.acceptance["missing_criteria"])
        self.assertIn("absence_bit_exports_equivariant_selector", self.acceptance["missing_criteria"])
        self.assertIn("new_inversion_odd_signed_source_exported", self.acceptance["missing_criteria"])
        self.assertIn("p2721_polarity_coupling_exported", self.acceptance["missing_criteria"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_recommendation_and_docs_are_updated(self):
        self.assertIn("information about absence of selector", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2697-P2748", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2748/S1698", MD.read_text(encoding="utf-8"))
        self.assertIn("P2748/S1698", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2748/S1698", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2748", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
