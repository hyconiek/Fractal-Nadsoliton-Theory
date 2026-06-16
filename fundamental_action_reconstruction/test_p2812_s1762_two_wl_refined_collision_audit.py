import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
JSON_PATH = ROOT / "generated" / "p2812_s1762_two_wl_refined_collision_audit.json"
MD_PATH = ROOT / "generated" / "p2812_s1762_two_wl_refined_collision_audit.md"


class P2812TwoWLRefinedCollisionAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(ROOT / "p2812_s1762_two_wl_refined_collision_audit.py")], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["two_wl_refinement_audit"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2812_TWO_WL_REFINED_COLLISION_AUDIT_REFINES_BUT_STILL_OBSTRUCTED_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2811"], "P2811_LOCAL_MOTIF_REFINED_SOURCE_CANDIDATE_AUDIT_REFINES_BUT_STILL_OBSTRUCTED_NO_CLOSURE")

    def test_exact_two_wl_refinement_counts(self):
        self.assertEqual(self.audit["decoded_graph_count"], 16828)
        self.assertEqual(self.audit["p2811_refined_class_count"], 16691)
        self.assertEqual(self.audit["p2811_refined_collision_graph_count"], 269)
        self.assertEqual(self.audit["two_wl_computed_graph_count"], 269)
        self.assertEqual(self.audit["two_wl_refined_class_count"], 16749)
        self.assertEqual(self.audit["two_wl_collision_class_count"], 79)
        self.assertEqual(self.audit["two_wl_collision_graph_count"], 158)
        self.assertEqual(self.audit["two_wl_max_class_size"], 2)
        self.assertEqual(self.audit["remaining_defect_canonical_minus_two_wl"], 79)
        self.assertEqual(self.audit["defect_reduction_vs_p2811"], 58)

    def test_acceptance_boundaries(self):
        self.assertTrue(self.acceptance["accepted_as_two_wl_refinement_audit"])
        self.assertFalse(self.acceptance["accepted_as_complete_canonical_source_carrier"])
        self.assertFalse(self.acceptance["accepted_as_strict_source_law_or_coupling"])
        self.assertFalse(self.payload["decision"]["negative_export_flags"]["role_bearing_ltotal_promoted"])
        self.assertFalse(self.payload["decision"]["negative_export_flags"]["toe_closure_exported"])

    def test_written_documents_reference_guardrail(self):
        self.assertIn("P2812/S1762", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2812/S1762", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2812/S1762", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2812", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
