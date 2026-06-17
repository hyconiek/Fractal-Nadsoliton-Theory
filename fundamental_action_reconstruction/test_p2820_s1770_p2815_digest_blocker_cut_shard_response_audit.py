import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2820_s1770_p2815_digest_blocker_cut_shard_response_audit.py"
JSON_PATH = ROOT / "generated" / "p2820_s1770_p2815_digest_blocker_cut_shard_response_audit.json"
MD_PATH = ROOT / "generated" / "p2820_s1770_p2815_digest_blocker_cut_shard_response_audit.md"
MANIFEST_PATH = ROOT / "generated" / "p2820_s1770_blocker_cut_shard_compact_manifest.json"


class P2820P2815DigestBlockerCutShardResponseAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2820_P2815_DIGEST_BLOCKER_CUT_SHARD_WITNESS_NO_FULL_COUPLING_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2818"], "P2818_LOCAL_EDGE_VARIATIONAL_RESPONSE_ENERGY_OBSTRUCTED_NO_COUPLING_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2819"], "P2819_P2815_DIGEST_BLOCKER_CUT_SAMPLE_WITNESS_NO_FULL_COUPLING_NO_CLOSURE")

    def test_shard_response_counts(self):
        audit = self.payload["p2815_digest_blocker_cut_shard_audit"]
        self.assertEqual(audit["decoded_graph_count"], 16828)
        self.assertEqual(audit["p2818_blocker_class_size"], 788)
        self.assertEqual(audit["shard_count"], 4)
        self.assertEqual(audit["shard_size"], 24)
        self.assertEqual(audit["computed_graph_count"], 96)
        self.assertEqual(audit["sample_fraction_of_blocker_cut"], [96, 788])
        self.assertEqual(audit["toggle_digest_refined_class_count"], 96)
        self.assertEqual(audit["toggle_digest_collision_class_count"], 0)
        self.assertEqual(audit["toggle_digest_collision_graph_count"], 0)
        self.assertEqual(audit["toggle_digest_max_class_size"], 1)
        self.assertEqual(audit["blocker_defect_after_toggle_digest"], 692)
        self.assertEqual(len(self.manifest["rows"]), 96)

    def test_acceptance_boundaries(self):
        matrix = self.payload["acceptance_matrix"]
        self.assertTrue(matrix["facts"]["p2818_local_response_rejected"])
        self.assertTrue(matrix["facts"]["p2819_sample_witness_positive"])
        self.assertTrue(matrix["facts"]["four_shards_audited"])
        self.assertTrue(matrix["facts"]["blocker_cut_is_fully_separated"])
        self.assertFalse(matrix["facts"]["entire_largest_blocker_cut_audited"])
        self.assertFalse(matrix["facts"]["all_16828_graphs_audited"])
        self.assertFalse(matrix["facts"]["all_p2818_collision_classes_audited"])
        self.assertFalse(matrix["facts"]["strict_graph_source_law_exported"])
        self.assertFalse(matrix["facts"]["typed_graph_to_K_or_Ltotal_coupling_exported"])
        self.assertTrue(matrix["accepted_as_four_shard_response_witness"])
        self.assertFalse(matrix["accepted_as_full_carrier_source_law_coupling"])
        self.assertTrue(matrix["accepted_as_bounded_no_closure_audit"])

    def test_documents_updated(self):
        self.assertIn("P2820/S1770", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2820/S1770", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2820/S1770", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2820", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
