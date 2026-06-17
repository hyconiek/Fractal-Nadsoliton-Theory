import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2830_s1780_p2815_digest_full_carrier_cached_audit.py"
JSON_PATH = ROOT / "generated" / "p2830_s1780_p2815_digest_full_carrier_cached_audit.json"
MD_PATH = ROOT / "generated" / "p2830_s1780_p2815_digest_full_carrier_cached_audit.md"
MANIFEST_PATH = ROOT / "generated" / "p2830_s1780_full_carrier_cached_digest_manifest.json"

SINGLETON_INDICES = [0, 1, 2, 3, 4, 5, 6, 7, 14, 15, 17, 18, 19, 20, 34, 61, 63, 66, 75, 82, 88, 89, 90, 103, 156, 161, 182, 186, 190, 192, 193, 197, 202, 205, 216, 230, 247, 255, 265, 272, 436, 438, 452, 478, 1513, 1673, 2050, 2051, 2277, 2289, 2308, 2309, 2375, 2452, 2571, 2625, 2638, 3557, 3560, 5133, 5525, 8771, 12774, 12878, 12879, 13463, 13536, 15533, 15536, 15542, 16807, 16827]


class P2830P2815DigestFullCarrierCachedAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2830_P2815_DIGEST_FULL_16828_CARRIER_SEPARATION_WITNESS_NO_SOURCE_LAW_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2818"], "P2818_LOCAL_EDGE_VARIATIONAL_RESPONSE_ENERGY_OBSTRUCTED_NO_COUPLING_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2829"], "P2829_P2815_DIGEST_RANKS_161_TO_272_BLOCKER_CUTS_WITNESS_NO_FULL_COUPLING_NO_CLOSURE")

    def test_full_carrier_counts(self):
        audit = self.payload["p2815_digest_full_carrier_cached_audit"]
        self.assertEqual(audit["decoded_graph_count"], 16828)
        self.assertEqual(audit["cached_collision_manifest_count"], 8)
        self.assertEqual(audit["cached_collision_row_count"], 16756)
        self.assertEqual(audit["cached_collision_unique_graph_count"], 16756)
        self.assertEqual(audit["p2818_collision_class_count_after_p2829"], 272)
        self.assertEqual(audit["p2818_collision_graph_count_after_p2829"], 16756)
        self.assertEqual(audit["fresh_singleton_p2818_local_response_class_count"], 72)
        self.assertEqual(audit["fresh_singleton_graph_count"], 72)
        self.assertEqual(audit["fresh_singleton_indices_0_based"], SINGLETON_INDICES)
        self.assertEqual(audit["full_carrier_row_count"], 16828)
        self.assertEqual(audit["full_carrier_unique_graph_count"], 16828)
        self.assertEqual(audit["full_carrier_missing_graph_count"], 0)
        self.assertEqual(audit["full_carrier_duplicate_graph_count"], 0)
        self.assertEqual(audit["full_carrier_digest_refined_class_count"], 16828)
        self.assertEqual(audit["full_carrier_digest_collision_class_count"], 0)
        self.assertEqual(audit["full_carrier_digest_collision_graph_count"], 0)
        self.assertEqual(audit["full_carrier_digest_max_class_size"], 1)
        self.assertEqual(audit["full_carrier_defect_after_digest"], 0)
        self.assertEqual(audit["full_carrier_class_size_histogram"], {"1": 16828})

    def test_manifest(self):
        self.assertEqual(self.manifest["cached_collision_row_count"], 16756)
        self.assertEqual(self.manifest["fresh_singleton_row_count"], 72)
        self.assertEqual(self.manifest["fresh_singleton_indices_0_based"], SINGLETON_INDICES)
        self.assertEqual(self.manifest["full_carrier_row_count"], 16828)
        self.assertEqual(len(self.manifest["cached_collision_manifests"]), 8)
        self.assertEqual(len(self.manifest["full_carrier_digest_sha256_by_graph_index"]), 16828)

    def test_acceptance_boundaries(self):
        matrix = self.payload["acceptance_matrix"]
        self.assertTrue(matrix["facts"]["p2818_local_response_rejected"])
        self.assertTrue(matrix["facts"]["p2829_all_collision_classes_witness_positive"])
        self.assertTrue(matrix["facts"]["all_p2818_collision_classes_audited"])
        self.assertTrue(matrix["facts"]["all_16828_graphs_audited"])
        self.assertTrue(matrix["facts"]["every_graph_index_present_once"])
        self.assertTrue(matrix["facts"]["full_carrier_digest_fully_separated"])
        self.assertFalse(matrix["facts"]["strict_graph_source_law_exported"])
        self.assertFalse(matrix["facts"]["typed_graph_to_K_or_Ltotal_coupling_exported"])
        self.assertFalse(matrix["facts"]["units_and_variational_derivative_exported"])
        self.assertTrue(matrix["accepted_as_full_16828_carrier_digest_separation_witness"])
        self.assertFalse(matrix["accepted_as_full_carrier_source_law_coupling"])
        self.assertTrue(matrix["accepted_as_bounded_no_closure_audit"])

    def test_documents_updated(self):
        self.assertIn("P2830/S1780", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2830/S1780", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2830/S1780", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2830", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
