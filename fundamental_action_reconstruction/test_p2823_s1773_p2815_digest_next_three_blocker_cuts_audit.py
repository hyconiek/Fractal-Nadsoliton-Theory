import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2823_s1773_p2815_digest_next_three_blocker_cuts_audit.py"
JSON_PATH = ROOT / "generated" / "p2823_s1773_p2815_digest_next_three_blocker_cuts_audit.json"
MD_PATH = ROOT / "generated" / "p2823_s1773_p2815_digest_next_three_blocker_cuts_audit.md"
MANIFEST_PATH = ROOT / "generated" / "p2823_s1773_next_three_blocker_cuts_compact_manifest.json"


class P2823P2815DigestNextThreeBlockerCutsAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2823_P2815_DIGEST_NEXT_THREE_BLOCKER_CUTS_WITNESS_NO_FULL_COUPLING_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2818"], "P2818_LOCAL_EDGE_VARIATIONAL_RESPONSE_ENERGY_OBSTRUCTED_NO_COUPLING_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2822"], "P2822_P2815_DIGEST_TOP_TWO_BLOCKER_CUTS_WITNESS_NO_FULL_COUPLING_NO_CLOSURE")

    def test_next_three_counts(self):
        audit = self.payload["p2815_digest_next_three_blocker_cuts_audit"]
        self.assertEqual(audit["decoded_graph_count"], 16828)
        self.assertEqual(audit["audited_rank_range_1_based"], [3, 5])
        self.assertEqual(audit["audited_p2818_collision_class_count"], 3)
        self.assertEqual(audit["audited_blocker_class_sizes"], [728, 691, 680])
        self.assertEqual(audit["audited_graph_count"], 2099)
        self.assertEqual(audit["previous_p2822_audited_collision_class_count"], 2)
        self.assertEqual(audit["previous_p2822_audited_graph_count"], 1564)
        self.assertEqual(audit["cumulative_audited_collision_class_count"], 5)
        self.assertEqual(audit["cumulative_audited_graph_count"], 3663)
        self.assertEqual(audit["combined_toggle_digest_refined_class_count"], 2099)
        self.assertEqual(audit["combined_toggle_digest_collision_class_count"], 0)
        self.assertEqual(audit["combined_toggle_digest_collision_graph_count"], 0)
        self.assertEqual(audit["combined_toggle_digest_max_class_size"], 1)
        self.assertEqual(audit["combined_defect_after_toggle_digest"], 0)
        self.assertEqual(len(self.manifest["rows"]), 2099)
        self.assertEqual(self.manifest["computed_new_row_count"], 2099)

    def test_per_class_counts(self):
        per_class = self.payload["p2815_digest_next_three_blocker_cuts_audit"]["per_class_results"]
        self.assertEqual([row["p2818_blocker_rank_1_based"] for row in per_class], [3, 4, 5])
        self.assertEqual([row["p2818_blocker_class_size"] for row in per_class], [728, 691, 680])
        self.assertEqual([row["toggle_digest_refined_class_count"] for row in per_class], [728, 691, 680])
        self.assertEqual([row["toggle_digest_collision_class_count"] for row in per_class], [0, 0, 0])
        self.assertEqual([row["blocker_defect_after_toggle_digest"] for row in per_class], [0, 0, 0])

    def test_acceptance_boundaries(self):
        matrix = self.payload["acceptance_matrix"]
        self.assertTrue(matrix["facts"]["p2818_local_response_rejected"])
        self.assertTrue(matrix["facts"]["p2822_top_two_witness_positive"])
        self.assertTrue(matrix["facts"]["expected_rank_range_audited"])
        self.assertTrue(matrix["facts"]["next_three_p2818_collision_classes_audited"])
        self.assertTrue(matrix["facts"]["each_audited_blocker_cut_fully_separated"])
        self.assertTrue(matrix["facts"]["combined_next_three_digest_fully_separated"])
        self.assertFalse(matrix["facts"]["all_p2818_collision_classes_audited"])
        self.assertFalse(matrix["facts"]["all_16828_graphs_audited"])
        self.assertFalse(matrix["facts"]["strict_graph_source_law_exported"])
        self.assertFalse(matrix["facts"]["typed_graph_to_K_or_Ltotal_coupling_exported"])
        self.assertTrue(matrix["accepted_as_next_three_blocker_cut_response_witness"])
        self.assertFalse(matrix["accepted_as_full_carrier_source_law_coupling"])
        self.assertTrue(matrix["accepted_as_bounded_no_closure_audit"])

    def test_documents_updated(self):
        self.assertIn("P2823/S1773", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2823/S1773", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2823/S1773", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2823", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
