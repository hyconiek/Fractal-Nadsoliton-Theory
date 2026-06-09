from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2630_s1580_strict_beta_source_vs_bridge_zbeta_separation.py"
OUT = ROOT / "generated" / "p2630_s1580_strict_beta_source_vs_bridge_zbeta_separation.json"
MD = ROOT / "generated" / "p2630_s1580_strict_beta_source_vs_bridge_zbeta_separation.md"


class P2630StrictBetaSourceVsBridgeZbetaSeparationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_is_present(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        self.assertTrue(all(data["count"] > 0 for data in audit["patterns"].values()))

    def test_truth_table_has_one_accepting_row(self) -> None:
        table = self.payload["bridge_source_truth_table"]
        self.assertEqual(table["row_count"], 8)
        self.assertEqual(table["accepting_row_count"], 1)
        self.assertTrue(table["unique_accepting_row"]["strict_internal_beta_source"])
        self.assertTrue(table["unique_accepting_row"]["legacy_uv_normalization_source"])
        self.assertTrue(table["unique_accepting_row"]["normalization_invariant_match_beta_micro_over_beta_strict_equals_1"])

    def test_current_assignment_rejects_even_if_p2603_granted(self) -> None:
        table = self.payload["bridge_source_truth_table"]
        current = table["current_assignment_after_p2629_even_if_p2603_granted"]
        self.assertTrue(current["strict_internal_beta_source"])
        self.assertFalse(current["legacy_uv_normalization_source"])
        self.assertFalse(current["normalization_invariant_match_beta_micro_over_beta_strict_equals_1"])
        self.assertFalse(table["current_accepts_bridge_positive_zbeta_source"])

    def test_source_type_matrix_blocks_reuse(self) -> None:
        matrix = self.payload["source_type_matrix"]
        self.assertEqual(len(matrix), 3)
        self.assertTrue(all(not row["can_supply_p2625_positive_zbeta_source"] for row in matrix))
        self.assertIn("P2603-style strict beta sourcehood", self.payload["bridge_source_truth_table"]["theorem"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))

    def test_recommendation_and_docs_updated(self) -> None:
        self.assertIn("new typed legacy/UV normalization theorem", self.payload["recommended_next_honest_step"])
        self.assertIn("P2630/S1580", MD.read_text(encoding="utf-8"))
        self.assertIn("P2630/S1580", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2630/S1580", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
