import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2827_s1777_p2815_digest_ranks_41_to_80_blocker_cuts_audit.py"
JSON_PATH = ROOT / "generated" / "p2827_s1777_p2815_digest_ranks_41_to_80_blocker_cuts_audit.json"
MD_PATH = ROOT / "generated" / "p2827_s1777_p2815_digest_ranks_41_to_80_blocker_cuts_audit.md"
MANIFEST_PATH = ROOT / "generated" / "p2827_s1777_ranks_41_to_80_blocker_cuts_compact_manifest.json"


class P2827P2815DigestRanks41To80BlockerCutsAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2827_P2815_DIGEST_RANKS_41_TO_80_BLOCKER_CUTS_WITNESS_NO_FULL_COUPLING_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2818"], "P2818_LOCAL_EDGE_VARIATIONAL_RESPONSE_ENERGY_OBSTRUCTED_NO_COUPLING_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2826"], "P2826_P2815_DIGEST_RANKS_21_TO_40_BLOCKER_CUTS_WITNESS_NO_FULL_COUPLING_NO_CLOSURE")

    def test_ranks_41_to_80_counts(self):
        audit = self.payload["p2815_digest_ranks_41_to_80_blocker_cuts_audit"]
        self.assertEqual(audit["decoded_graph_count"], 16828)
        self.assertEqual(audit["audited_rank_range_1_based"], [41, 80])
        self.assertEqual(audit["audited_p2818_collision_class_count"], 40)
        self.assertEqual(audit["audited_blocker_class_sizes"], [104, 96, 93, 85, 82, 80, 79, 79, 78, 75, 69, 67, 66, 66, 63, 61, 60, 60, 58, 57, 55, 54, 53, 52, 52, 45, 42, 42, 41, 40, 37, 37, 37, 36, 34, 34, 33, 29, 28, 28])
        self.assertEqual(audit["audited_graph_count"], 2287)
        self.assertEqual(audit["previous_p2826_cumulative_audited_collision_class_count"], 40)
        self.assertEqual(audit["previous_p2826_cumulative_audited_graph_count"], 13136)
        self.assertEqual(audit["cumulative_audited_collision_class_count"], 80)
        self.assertEqual(audit["cumulative_audited_graph_count"], 15423)
        self.assertEqual(audit["combined_toggle_digest_refined_class_count"], 2287)
        self.assertEqual(audit["combined_toggle_digest_collision_class_count"], 0)
        self.assertEqual(audit["combined_toggle_digest_collision_graph_count"], 0)
        self.assertEqual(audit["combined_toggle_digest_max_class_size"], 1)
        self.assertEqual(audit["combined_defect_after_toggle_digest"], 0)
        self.assertEqual(len(self.manifest["rows"]), 2287)
        self.assertEqual(self.manifest["computed_new_row_count"], 2287)

    def test_per_class_counts(self):
        per_class = self.payload["p2815_digest_ranks_41_to_80_blocker_cuts_audit"]["per_class_results"]
        sizes = [104, 96, 93, 85, 82, 80, 79, 79, 78, 75, 69, 67, 66, 66, 63, 61, 60, 60, 58, 57, 55, 54, 53, 52, 52, 45, 42, 42, 41, 40, 37, 37, 37, 36, 34, 34, 33, 29, 28, 28]
        self.assertEqual([row["p2818_blocker_rank_1_based"] for row in per_class], list(range(41, 81)))
        self.assertEqual([row["p2818_blocker_class_size"] for row in per_class], sizes)
        self.assertEqual([row["toggle_digest_refined_class_count"] for row in per_class], sizes)
        self.assertEqual([row["toggle_digest_collision_class_count"] for row in per_class], [0] * 40)
        self.assertEqual([row["blocker_defect_after_toggle_digest"] for row in per_class], [0] * 40)

    def test_acceptance_boundaries(self):
        matrix = self.payload["acceptance_matrix"]
        self.assertTrue(matrix["facts"]["p2818_local_response_rejected"])
        self.assertTrue(matrix["facts"]["p2826_ranks_21_to_40_witness_positive"])
        self.assertTrue(matrix["facts"]["expected_rank_range_audited"])
        self.assertTrue(matrix["facts"]["ranks_41_to_80_p2818_collision_classes_audited"])
        self.assertTrue(matrix["facts"]["each_audited_blocker_cut_fully_separated"])
        self.assertTrue(matrix["facts"]["combined_ranks_41_to_80_digest_fully_separated"])
        self.assertFalse(matrix["facts"]["all_p2818_collision_classes_audited"])
        self.assertFalse(matrix["facts"]["all_16828_graphs_audited"])
        self.assertFalse(matrix["facts"]["strict_graph_source_law_exported"])
        self.assertFalse(matrix["facts"]["typed_graph_to_K_or_Ltotal_coupling_exported"])
        self.assertTrue(matrix["accepted_as_ranks_41_to_80_blocker_cut_response_witness"])
        self.assertFalse(matrix["accepted_as_full_carrier_source_law_coupling"])
        self.assertTrue(matrix["accepted_as_bounded_no_closure_audit"])

    def test_documents_updated(self):
        self.assertIn("P2827/S1777", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2827/S1777", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2827/S1777", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2827", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
