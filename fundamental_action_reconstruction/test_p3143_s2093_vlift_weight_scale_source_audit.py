import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3143_s2093_vlift_weight_scale_source_audit.py"
OUT = ROOT / "generated" / "p3143_s2093_vlift_weight_scale_source_audit.json"
MD = ROOT / "generated" / "p3143_s2093_vlift_weight_scale_source_audit.md"


class P3143VliftWeightScaleSourceAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_inputs_and_counts(self):
        self.assertEqual(self.payload["status"], "P3143_VLIFT_WEIGHT_SCALE_SOURCE_AUDIT_BOUNDED_NO_GO")
        self.assertTrue(self.payload["input_hashes"]["P3142"])
        counts = self.payload["finite_theorem"]["finite_counts"]
        self.assertEqual(counts["candidates_tested"], 5)
        self.assertEqual(counts["targets_tested"], 4)
        self.assertEqual(counts["candidate_target_rows"], 20)
        self.assertEqual(counts["positive_value_rows"], 16)
        self.assertEqual(counts["accepted_strict_source_rows"], 0)

    def test_obligation_matrix_blocks_all_rows(self):
        rows = self.payload["candidate_target_rows"]
        self.assertEqual(len(rows), 20)
        self.assertTrue(all(not row["accepted_as_strict_source"] for row in rows))
        strict_source_defects = [row for row in rows if not row["obligations"]["strict_source_law"]]
        unit_defects = [row for row in rows if not row["obligations"]["unit_bearing"]]
        self.assertEqual(len(strict_source_defects), 20)
        self.assertEqual(len(unit_defects), 20)

    def test_defect_table_and_recommendation(self):
        defects = {row["obligation"]: row for row in self.payload["source_defect_table"]}
        self.assertEqual(defects["positive scalar value"]["rows_passing"], 16)
        self.assertEqual(defects["strict source law"]["rows_passing"], 0)
        self.assertEqual(defects["unit-bearing normalization"]["rows_passing"], 0)
        self.assertIn("Upsilon_sel", self.payload["decision"]["next_honest_step"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_docs_updated(self):
        self.assertIn("P3143/S2093", MD.read_text(encoding="utf-8"))
        self.assertIn("P3143/S2093", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3143/S2093", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3143", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
