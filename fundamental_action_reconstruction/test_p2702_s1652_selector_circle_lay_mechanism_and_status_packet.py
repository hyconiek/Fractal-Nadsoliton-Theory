from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2702_s1652_selector_circle_lay_mechanism_and_status_packet.py"
OUT = ROOT / "generated" / "p2702_s1652_selector_circle_lay_mechanism_and_status_packet.json"
MD = ROOT / "generated" / "p2702_s1652_selector_circle_lay_mechanism_and_status_packet.md"


class P2702SelectorCircleLayMechanismStatusTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_grep_covers_selector_circle_language(self) -> None:
        audit = self.payload["content_grep"]
        self.assertIn("P2702", audit["mode"])
        for key in ("selector_status_chain", "circle_language", "premise_vs_strict", "convexity_language"):
            self.assertIn(key, audit["patterns"])
            self.assertGreater(audit["patterns"][key]["count"], 0)

    def test_circle_models_separate_premise_selector_from_strict_candidate(self) -> None:
        circle = self.payload["z12_circle_models"]
        self.assertEqual(circle["units"], [1, 5, 7, 11])
        self.assertEqual(circle["generator_orbit_of_plus_one"], [1, 5, 7, 11])
        self.assertTrue(circle["minus_one_in_plus_one_orbit"])
        self.assertTrue(circle["premise_selector_breaks_symmetry_in_toy_model"])
        self.assertFalse(circle["aut_invariant_candidate_breaks_symmetry_in_toy_model"])

    def test_state_reads_and_lay_rows_preserve_no_false_pass(self) -> None:
        state = self.payload["state_reads"]
        self.assertTrue(state["p2698_real_direction_support_preserved"])
        self.assertTrue(state["p2699_aut_invariant_no_go"])
        self.assertTrue(state["p2700_exhaustive_no_go"])
        self.assertTrue(state["p2701_no_provider"])
        rows = self.payload["lay_explanation_rows"]
        self.assertEqual(len(rows), 4)
        self.assertIn("selector", rows[0]["plain_answer"])
        self.assertIn("convexity", rows[2]["question"].lower())

    def test_decision_and_docs(self) -> None:
        decision = self.payload["decision"]
        self.assertTrue(decision["premise_selector_can_pick_one_direction_in_toy_model"])
        self.assertFalse(decision["strict_aut_invariant_candidate_picks_one_direction"])
        self.assertTrue(decision["no_new_closure_exported"])
        self.assertIn("no-new-live-frontier", decision["next_honest_step"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2702/S1652", MD.read_text(encoding="utf-8"))
        self.assertIn("P2702/S1652", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2702/S1652", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2702/S1652", (ROOT.parent / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
