from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2700_s1650_exhaustive_aut_invariant_selector_functional_enumeration_no_go.py"
OUT = ROOT / "generated" / "p2700_s1650_exhaustive_aut_invariant_selector_functional_enumeration_no_go.json"
MD = ROOT / "generated" / "p2700_s1650_exhaustive_aut_invariant_selector_functional_enumeration_no_go.md"


class P2700ExhaustiveAutInvariantSelectorFunctionalEnumerationNoGoTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_grep_covers_exhaustive_selector_functional_boundary(self) -> None:
        audit = self.payload["content_grep"]
        self.assertIn("P2700", audit["mode"])
        for key in ("p2699_no_go", "selector_functional_terms", "z12_aut_terms", "qw2191_boundaries", "forbidden_promotions"):
            self.assertIn(key, audit["patterns"])
            self.assertGreater(audit["patterns"][key]["count"], 0)

    def test_state_reads_preserve_p2699_and_nonexport_boundaries(self) -> None:
        state = self.payload["state_reads"]
        self.assertTrue(state["p2699_bounded_no_go"])
        self.assertTrue(state["p2699_requested_new_nonpremise_source"])
        self.assertTrue(state["h39_qw2191_still_open"])
        self.assertTrue(state["p739_pair12_upgrade_unexported"])
        self.assertTrue(state["p740_pair12_upgrade_unexported"])

    def test_exhaustive_calculation_counts_and_no_selector(self) -> None:
        calc = self.payload["exhaustive_calculation"]
        self.assertEqual(calc["units"], [1, 5, 7, 11])
        self.assertEqual(len(calc["orbits"]), 6)
        self.assertTrue(calc["complete_boolean_space_checked"])
        self.assertEqual(calc["boolean_predicate_enumeration"]["predicate_count"], 64)
        self.assertEqual(calc["boolean_predicate_enumeration"]["unique_unit_singleton_hit_count"], 0)
        self.assertEqual(calc["boolean_predicate_enumeration"]["plus_without_minus_hit_count"], 0)
        self.assertEqual(calc["orbit_constant_score_enumeration"]["score_function_count"], 15625)
        self.assertEqual(calc["orbit_constant_score_enumeration"]["unique_unit_argmax_hit_count"], 0)
        self.assertEqual(calc["orbit_constant_score_enumeration"]["plus_minus_separating_score_count"], 0)
        self.assertFalse(calc["selector_found"])
        self.assertFalse(calc["plus_minus_distinction_found"])

    def test_obstruction_rows_and_decision_preserve_no_false_pass(self) -> None:
        rows = self.payload["obstruction_rows"]
        self.assertEqual(len(rows), 5)
        self.assertTrue(all(not row["passes"] for row in rows))
        decision = self.payload["decision"]
        self.assertTrue(decision["bounded_no_go_now"])
        self.assertIn("no-new-live-frontier", decision["next_honest_step"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2700/S1650", MD.read_text(encoding="utf-8"))
        self.assertIn("P2700/S1650", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2700/S1650", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2700/S1650", (ROOT.parent / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
