import csv
import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
JSON_PATH = ROOT / "generated" / "p2806_s1756_girth4_reproducible_record_label_dataset_audit.json"
MD_PATH = ROOT / "generated" / "p2806_s1756_girth4_reproducible_record_label_dataset_audit.md"
CSV_PATH = ROOT / "generated" / "p2806_s1756_girth4_record_label_dataset.csv"


class P2806Girth4ReproducibleRecordLabelDatasetAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(ROOT / "p2806_s1756_girth4_reproducible_record_label_dataset_audit.py")], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["record_label_dataset_audit"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_status_and_p2805_input(self):
        self.assertEqual(
            self.payload["status"],
            "P2806_GIRTH4_REPRODUCIBLE_RECORD_LABEL_DATASET_UNIQUE_RECORDS_NO_ISOMORPHISM_CANONICAL_LABELS_NO_SOURCE_LAW_NO_CLOSURE",
        )
        self.assertEqual(
            self.payload["input_statuses"]["P2805"],
            "P2805_GIRTH4_SPECTRAL_COLLISION_ISOMORPHISM_REFINEMENT_NO_CANONICAL_LABELS_NO_SOURCE_LAW_NO_CLOSURE",
        )
        self.assertEqual(self.audit["p2805_resolved_total_isomorphism_classes_after_refinement"], 16828)

    def test_record_label_dataset_counts(self):
        self.assertEqual(self.audit["decoded_graph_count"], 16828)
        self.assertEqual(self.audit["record_label_count"], 16828)
        self.assertEqual(self.audit["unique_record_graph6_label_count"], 16828)
        self.assertEqual(self.audit["unique_record_graph6_sha256_count"], 16828)
        self.assertEqual(self.audit["duplicate_record_label_count"], 0)
        self.assertEqual(self.audit["duplicate_record_hash_count"], 0)
        self.assertTrue(CSV_PATH.exists())

    def test_csv_shape_and_hash_samples(self):
        with CSV_PATH.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 16828)
        self.assertEqual(rows[0]["graph_index_1_based"], "1")
        self.assertEqual(rows[-1]["graph_index_1_based"], "16828")
        self.assertEqual(rows[0]["record_graph6_label"], self.audit["sample_rows"][0]["record_graph6_label"])
        self.assertEqual(rows[-1]["record_graph6_sha256"], self.audit["tail_rows"][-1]["record_graph6_sha256"])

    def test_acceptance_boundaries(self):
        self.assertTrue(self.acceptance["accepted_as_reproducible_record_label_dataset"])
        self.assertFalse(self.acceptance["accepted_as_isomorphism_canonical_label_dataset"])
        self.assertFalse(self.acceptance["accepted_as_independent_canonical_label_crosscheck"])
        self.assertFalse(self.acceptance["accepted_as_strict_spectral_source_law"])
        self.assertFalse(self.acceptance["accepted_as_ltotal_or_toe_promotion"])
        self.assertIn("isomorphism_canonical_label_dataset_exported", self.acceptance["missing_criteria"])
        self.assertFalse(self.payload["decision"]["negative_export_flags"]["toe_closure_exported"])

    def test_written_documents_reference_guardrail(self):
        self.assertIn("P2806/S1756", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2806/S1756", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2806/S1756", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2806", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
