import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2816_s1766_edge_toggle_source_law_coupling_acceptance_audit.py"
JSON_PATH = ROOT / "generated" / "p2816_s1766_edge_toggle_source_law_coupling_acceptance_audit.json"
MD_PATH = ROOT / "generated" / "p2816_s1766_edge_toggle_source_law_coupling_acceptance_audit.md"


class P2816EdgeToggleSourceLawCouplingAcceptanceAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2816_EDGE_TOGGLE_SOURCE_FUNCTIONAL_REJECTED_NO_COUPLING_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2815"], "P2815_EDGE_TOGGLE_RESPONSE_SEPARATES_CARRIER_NO_SOURCE_LAW_NO_CLOSURE")

    def test_rank_candidate_counts(self):
        audit = self.payload["q_rank_source_candidate_audit"]
        self.assertEqual(audit["decoded_graph_count"], 16828)
        self.assertEqual(audit["input_p2815_refined_class_count"], 16828)
        self.assertEqual(audit["input_p2815_collision_class_count"], 0)
        self.assertEqual(audit["q_rank_class_count"], 16828)
        self.assertEqual(audit["q_rank_collision_class_count"], 0)
        self.assertEqual(audit["q_rank_collision_graph_count"], 0)
        self.assertEqual(audit["remaining_defect_canonical_minus_q_rank"], 0)
        self.assertEqual(audit["rank_domain_min"], 0)
        self.assertEqual(audit["rank_domain_max"], 16827)

    def test_acceptance_boundaries(self):
        matrix = self.payload["acceptance_matrix"]
        self.assertTrue(matrix["facts"]["p2815_complete_carrier_available"])
        self.assertTrue(matrix["facts"]["functional_separates_full_carrier"])
        self.assertFalse(matrix["facts"]["dimensionful_normalization_exported"])
        self.assertFalse(matrix["facts"]["strict_graph_source_law_exported"])
        self.assertFalse(matrix["facts"]["typed_graph_to_K_or_Ltotal_coupling_exported"])
        self.assertFalse(matrix["accepted_as_source_law_coupling"])
        self.assertTrue(matrix["accepted_as_bounded_candidate_rejection"])

    def test_documents_updated(self):
        self.assertIn("P2816/S1766", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2816/S1766", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2816/S1766", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2816", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
