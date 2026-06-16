import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
JSON_PATH = ROOT / "generated" / "p2814_s1764_exact_automorphism_orbit_residual_audit.json"
MD_PATH = ROOT / "generated" / "p2814_s1764_exact_automorphism_orbit_residual_audit.md"


class P2814ExactAutomorphismOrbitResidualAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(ROOT / "p2814_s1764_exact_automorphism_orbit_residual_audit.py")], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["automorphism_orbit_audit"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2814_EXACT_AUTOMORPHISM_ORBIT_RESIDUAL_OBSTRUCTION_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2813"], "P2813_THREE_WL_HISTOGRAM_RESIDUAL_COLLISION_OBSTRUCTION_NO_CLOSURE")

    def test_exact_automorphism_orbit_counts(self):
        self.assertEqual(self.audit["decoded_graph_count"], 16828)
        self.assertEqual(self.audit["p2813_refined_class_count"], 16749)
        self.assertEqual(self.audit["p2813_collision_class_count"], 79)
        self.assertEqual(self.audit["p2813_collision_graph_count"], 158)
        self.assertEqual(self.audit["automorphism_orbit_computed_graph_count"], 158)
        self.assertEqual(self.audit["automorphism_orbit_refined_class_count"], 16771)
        self.assertEqual(self.audit["automorphism_orbit_collision_class_count"], 57)
        self.assertEqual(self.audit["automorphism_orbit_collision_graph_count"], 114)
        self.assertEqual(self.audit["automorphism_orbit_max_class_size"], 2)
        self.assertEqual(self.audit["remaining_defect_canonical_minus_automorphism_orbit"], 57)
        self.assertEqual(self.audit["defect_reduction_vs_p2813"], 22)
        self.assertEqual(self.audit["truncated_automorphism_search_count"], 0)
        self.assertEqual(self.audit["automorphism_group_order_histogram_on_residual"], {"1": 112, "2": 40, "4": 6})

    def test_acceptance_boundaries(self):
        self.assertTrue(self.acceptance["accepted_as_exact_automorphism_orbit_obstruction_audit"])
        self.assertFalse(self.acceptance["accepted_as_complete_canonical_source_carrier"])
        self.assertFalse(self.acceptance["accepted_as_strict_source_law_or_coupling"])
        self.assertFalse(self.payload["decision"]["negative_export_flags"]["role_bearing_ltotal_promoted"])
        self.assertFalse(self.payload["decision"]["negative_export_flags"]["toe_closure_exported"])

    def test_written_documents_reference_guardrail(self):
        self.assertIn("P2814/S1764", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2814/S1764", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2814/S1764", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2814", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
