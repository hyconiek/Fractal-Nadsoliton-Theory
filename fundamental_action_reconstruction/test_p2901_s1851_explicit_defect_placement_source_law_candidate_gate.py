import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2901_s1851_explicit_defect_placement_source_law_candidate_gate.py"
JSON_PATH = ROOT / "generated" / "p2901_s1851_explicit_defect_placement_source_law_candidate_gate.json"
MD_PATH = ROOT / "generated" / "p2901_s1851_explicit_defect_placement_source_law_candidate_gate.md"


class P2901ExplicitDefectPlacementSourceLawCandidateGateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.acceptance = cls.payload["acceptance_matrix"]
        cls.constructed = cls.payload["constructed_theoretical_object"]
        cls.flags = cls.payload["decision"]["negative_export_flags"]

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2901_EXPLICIT_DEFECT_PLACEMENT_SOURCE_LAW_CANDIDATE_CONDITIONAL_NO_CLOSURE")
        self.assertTrue(self.acceptance["p2900_rechecked"])
        self.assertTrue(self.acceptance["formal_candidate_constructed"])
        self.assertTrue(self.acceptance["has_correct_formal_defect_and_9_over_5_shape"])

    def test_candidate_family_remains_imported(self):
        self.assertEqual(self.acceptance["candidate_parameter_count"], 24)
        self.assertEqual(self.acceptance["translation_orbit_count"], 2)
        self.assertEqual(self.acceptance["translation_orbit_sizes"], [12, 12])
        self.assertEqual(self.acceptance["unique_unimported_choice_count"], 0)
        self.assertFalse(self.acceptance["nonimported_basepoint_supplied_by_formula"])
        self.assertFalse(self.acceptance["nonimported_polarity_supplied_by_formula"])
        self.assertFalse(self.acceptance["accepted_as_missing_object"])

    def test_constructed_edges_have_9_over_5_offsets(self):
        rows = self.constructed["constructed_rows"]
        self.assertEqual(len(rows), 24)
        offsets = {row["density_support"]["carrier_offset_mod_12"] for row in rows}
        self.assertEqual(offsets, {5, 7})
        self.assertTrue(all(row["signal_zero_at_basepoint"] for row in rows))

    def test_false_export_flags(self):
        self.assertFalse(any(self.flags.values()))
        self.assertFalse(self.flags["toe_closure_exported"])
        self.assertFalse(self.flags["ltotal_exported"])

    def test_documents_updated(self):
        self.assertIn("P2901/S1851", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2901/S1851", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2901/S1851", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2901", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
