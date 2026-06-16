import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
JSON_PATH = ROOT / "generated" / "p2815_s1765_edge_toggle_response_residual_audit.json"
MD_PATH = ROOT / "generated" / "p2815_s1765_edge_toggle_response_residual_audit.md"


class P2815EdgeToggleResponseResidualAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(ROOT / "p2815_s1765_edge_toggle_response_residual_audit.py")], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["edge_toggle_response_audit"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2815_EDGE_TOGGLE_RESPONSE_SEPARATES_CARRIER_NO_SOURCE_LAW_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2814"], "P2814_EXACT_AUTOMORPHISM_ORBIT_RESIDUAL_OBSTRUCTION_NO_CLOSURE")

    def test_exact_edge_toggle_counts(self):
        self.assertEqual(self.audit["decoded_graph_count"], 16828)
        self.assertEqual(self.audit["p2814_refined_class_count"], 16771)
        self.assertEqual(self.audit["p2814_collision_class_count"], 57)
        self.assertEqual(self.audit["p2814_collision_graph_count"], 114)
        self.assertEqual(self.audit["edge_toggle_computed_graph_count"], 114)
        self.assertEqual(self.audit["edge_toggle_refined_class_count"], 16828)
        self.assertEqual(self.audit["edge_toggle_collision_class_count"], 0)
        self.assertEqual(self.audit["edge_toggle_collision_graph_count"], 0)
        self.assertEqual(self.audit["edge_toggle_max_class_size"], 1)
        self.assertEqual(self.audit["remaining_defect_canonical_minus_edge_toggle"], 0)
        self.assertEqual(self.audit["defect_reduction_vs_p2814"], 57)
        self.assertEqual(self.audit["computed_signature_distinct_count"], 114)

    def test_acceptance_boundaries(self):
        self.assertTrue(self.acceptance["accepted_as_edge_toggle_response_audit"])
        self.assertTrue(self.acceptance["accepted_as_complete_canonical_source_carrier"])
        self.assertFalse(self.acceptance["accepted_as_strict_source_law_or_coupling"])
        self.assertFalse(self.payload["decision"]["negative_export_flags"]["role_bearing_ltotal_promoted"])
        self.assertFalse(self.payload["decision"]["negative_export_flags"]["toe_closure_exported"])

    def test_written_documents_reference_guardrail(self):
        self.assertIn("P2815/S1765", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2815/S1765", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2815/S1765", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2815", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
