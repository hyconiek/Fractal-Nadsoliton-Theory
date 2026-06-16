import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
JSON_PATH = ROOT / "generated" / "p2813_s1763_three_wl_refined_collision_audit.json"
MD_PATH = ROOT / "generated" / "p2813_s1763_three_wl_refined_collision_audit.md"


class P2813ThreeWLRefinedCollisionAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(ROOT / "p2813_s1763_three_wl_refined_collision_audit.py")], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["three_wl_refinement_audit"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2813_THREE_WL_HISTOGRAM_RESIDUAL_COLLISION_OBSTRUCTION_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2812"], "P2812_TWO_WL_REFINED_COLLISION_AUDIT_REFINES_BUT_STILL_OBSTRUCTED_NO_CLOSURE")

    def test_exact_three_wl_obstruction_counts(self):
        self.assertEqual(self.audit["decoded_graph_count"], 16828)
        self.assertEqual(self.audit["p2812_two_wl_refined_class_count"], 16749)
        self.assertEqual(self.audit["p2812_collision_class_count"], 79)
        self.assertEqual(self.audit["p2812_collision_graph_count"], 158)
        self.assertEqual(self.audit["three_wl_computed_graph_count"], 158)
        self.assertEqual(self.audit["three_wl_refined_class_count"], 16749)
        self.assertEqual(self.audit["three_wl_collision_class_count"], 79)
        self.assertEqual(self.audit["three_wl_collision_graph_count"], 158)
        self.assertEqual(self.audit["three_wl_max_class_size"], 2)
        self.assertEqual(self.audit["remaining_defect_canonical_minus_three_wl"], 79)
        self.assertEqual(self.audit["defect_reduction_vs_p2812"], 0)

    def test_acceptance_boundaries(self):
        self.assertTrue(self.acceptance["accepted_as_three_wl_residual_obstruction_audit"])
        self.assertFalse(self.acceptance["accepted_as_complete_canonical_source_carrier"])
        self.assertFalse(self.acceptance["accepted_as_strict_source_law_or_coupling"])
        self.assertFalse(self.payload["decision"]["negative_export_flags"]["role_bearing_ltotal_promoted"])
        self.assertFalse(self.payload["decision"]["negative_export_flags"]["toe_closure_exported"])

    def test_written_documents_reference_guardrail(self):
        self.assertIn("P2813/S1763", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2813/S1763", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2813/S1763", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2813", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
