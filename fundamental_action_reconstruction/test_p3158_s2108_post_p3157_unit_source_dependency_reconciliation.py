import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3158_s2108_post_p3157_unit_source_dependency_reconciliation.py"
OUT = ROOT / "generated" / "p3158_s2108_post_p3157_unit_source_dependency_reconciliation.json"
MD = ROOT / "generated" / "p3158_s2108_post_p3157_unit_source_dependency_reconciliation.md"


class P3158PostP3157UnitSourceDependencyReconciliationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_dependency_counts_and_status(self):
        self.assertEqual(self.payload["status"], "P3158_POST_P3157_UNIT_SOURCE_DEPENDENCY_RECONCILIATION_NO_STRICT_UNIT")
        counts = self.payload["finite_theorem"]["finite_counts"]
        self.assertEqual(counts["dag_nodes"], 14)
        self.assertGreater(counts["dag_edges"], 0)
        self.assertEqual(counts["exported_nodes"], 0)
        self.assertEqual(counts["closed_nodes_now"], 0)
        self.assertGreater(counts["missing_leaf_cut_size"], 0)

    def test_missing_leaf_cut_and_no_exports(self):
        self.assertIn("Lambda_origin_source_localizer", self.payload["missing_leaf_cut"])
        self.assertIn("positive_torsor_source_law", self.payload["missing_leaf_cut"])
        self.assertIn("Lambda_origin_source_localizer", self.payload["decision"]["next_honest_step"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))
        self.assertGreater(self.payload["content_grep"]["omega_dim_kdim"]["count"], 0)

    def test_docs_updated(self):
        self.assertIn("P3158/S2108", MD.read_text(encoding="utf-8"))
        self.assertIn("P3158/S2108", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3158/S2108", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3158", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
