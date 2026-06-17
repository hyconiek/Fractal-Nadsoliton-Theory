import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2817_s1767_structural_observable_source_law_obstruction_audit.py"
JSON_PATH = ROOT / "generated" / "p2817_s1767_structural_observable_source_law_obstruction_audit.json"
MD_PATH = ROOT / "generated" / "p2817_s1767_structural_observable_source_law_obstruction_audit.md"


class P2817StructuralObservableSourceLawObstructionAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2817_STRUCTURAL_OBSERVABLE_SOURCE_CANDIDATE_OBSTRUCTED_NO_COUPLING_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2816"], "P2816_EDGE_TOGGLE_SOURCE_FUNCTIONAL_REJECTED_NO_COUPLING_NO_CLOSURE")

    def test_structural_observable_counts(self):
        audit = self.payload["q_struct_source_candidate_audit"]
        self.assertEqual(audit["decoded_graph_count"], 16828)
        self.assertEqual(audit["q_struct_class_count"], 228)
        self.assertEqual(audit["q_struct_collision_class_count"], 165)
        self.assertEqual(audit["q_struct_collision_graph_count"], 16765)
        self.assertEqual(audit["q_struct_max_class_size"], 788)
        self.assertEqual(audit["remaining_defect_canonical_minus_q_struct"], 16600)

    def test_acceptance_boundaries(self):
        matrix = self.payload["acceptance_matrix"]
        self.assertTrue(matrix["facts"]["p2816_completed_rank_candidate_rejected"])
        self.assertTrue(matrix["facts"]["exactly_one_nonordinal_observable_tested"])
        self.assertTrue(matrix["facts"]["observable_has_explicit_dimensionless_normalization_candidate"])
        self.assertFalse(matrix["facts"]["observable_separates_full_carrier"])
        self.assertFalse(matrix["facts"]["strict_graph_source_law_exported"])
        self.assertFalse(matrix["facts"]["typed_graph_to_K_or_Ltotal_coupling_exported"])
        self.assertTrue(matrix["accepted_as_structural_observable_obstruction_audit"])
        self.assertFalse(matrix["accepted_as_source_law_coupling"])
        self.assertTrue(matrix["accepted_as_bounded_candidate_rejection"])

    def test_documents_updated(self):
        self.assertIn("P2817/S1767", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2817/S1767", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2817/S1767", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2817", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
