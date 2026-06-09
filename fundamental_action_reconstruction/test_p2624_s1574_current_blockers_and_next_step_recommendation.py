from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2624_s1574_current_blockers_and_next_step_recommendation.py"
OUT = ROOT / "generated" / "p2624_s1574_current_blockers_and_next_step_recommendation.json"
MD = ROOT / "generated" / "p2624_s1574_current_blockers_and_next_step_recommendation.md"


class P2624CurrentBlockersAndNextStepRecommendationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_is_present(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        self.assertTrue(all(data["count"] > 0 for data in audit["patterns"].values()))

    def test_blocker_truth_table_has_single_toe_accepting_row(self) -> None:
        blockers = self.payload["current_blocker_vector"]
        self.assertEqual(blockers["truth_table_row_count"], 32)
        self.assertEqual(blockers["toe_accepting_rows"], 1)
        self.assertEqual(blockers["ltotal_accepting_rows"], 4)
        self.assertFalse(blockers["current_row"]["strict_kernel_is_full_toe_kernel"])
        self.assertEqual(set(blockers["current_row"]["missing_for_toe"]), set(blockers["atoms"]))

    def test_strict_kernel_symptoms_are_not_full_kernel(self) -> None:
        answer = self.payload["answer_to_full_kernel_question"]
        self.assertIn("No", answer["short_answer"])
        self.assertIn("not an admitted full ToE kernel", answer["precise_status"])
        self.assertTrue(all(row["status"] != "full_kernel" for row in self.payload["strict_kernel_toe_symptom_matrix"]))

    def test_recommendation_targets_nonlinear_damping_not_selector_loop(self) -> None:
        rec = self.payload["recommended_next_honest_step"]
        self.assertIn("nonlinear damping", rec["recommended_next_packet"])
        self.assertIn("stop cycling", rec["recommended_focus"])
        self.assertIn("nonlinear_damping_completion_source", rec["recommended_focus"])
        self.assertTrue(any("do not call K_strict_gate" in item for item in rec["forbidden_shortcuts"]))

    def test_no_negative_exports_and_docs_updated(self) -> None:
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2624/S1574", MD.read_text(encoding="utf-8"))
        self.assertIn("P2624/S1574", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2624/S1574", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
