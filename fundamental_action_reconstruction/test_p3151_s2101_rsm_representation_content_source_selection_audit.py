import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3151_s2101_rsm_representation_content_source_selection_audit.py"
OUT = ROOT / "generated" / "p3151_s2101_rsm_representation_content_source_selection_audit.json"
MD = ROOT / "generated" / "p3151_s2101_rsm_representation_content_source_selection_audit.md"


class P3151RsmRepresentationContentSourceSelectionAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_nonunique_counts(self):
        self.assertEqual(self.payload["status"], "P3151_RSM_REPRESENTATION_CONTENT_SOURCE_SELECTION_MULTI_WITNESS_NO_STRICT_SOURCE")
        scan = self.payload["finite_scan"]
        self.assertEqual(scan["total_shapes_scanned"], 6 ** 6)
        self.assertTrue(scan["sm_shape_present"])
        self.assertFalse(scan["unique_shape_selected"])
        self.assertGreater(scan["coarse_anomaly_pass_count"], 1)
        self.assertGreater(scan["distinct_dimension_pattern_count"], 1)

    def test_no_strict_source_exported(self):
        counts = self.payload["finite_theorem"]["finite_counts"]
        self.assertEqual(counts["strict_accepted_source_rows"], 0)
        self.assertFalse(any(row["selects_sm_shape"] and row["strict_source_law"] and row["noncircular"] for row in self.payload["candidate_source_rows"]))
        self.assertIn("P3152", self.payload["decision"]["next_honest_step"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_docs_updated(self):
        self.assertIn("P3151/S2101", MD.read_text(encoding="utf-8"))
        self.assertIn("P3151/S2101", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3151/S2101", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3151", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
