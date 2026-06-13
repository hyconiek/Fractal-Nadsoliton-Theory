from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2699_s1649_z12_fractal_information_aut_invariant_selector_source_no_go.py"
OUT = ROOT / "generated" / "p2699_s1649_z12_fractal_information_aut_invariant_selector_source_no_go.json"
MD = ROOT / "generated" / "p2699_s1649_z12_fractal_information_aut_invariant_selector_source_no_go.md"


class P2699Z12FractalInformationAutInvariantSelectorSourceNoGoTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_grep_covers_p2698_followup_and_boundaries(self) -> None:
        audit = self.payload["content_grep"]
        self.assertIn("Aut(Z12)", audit["mode"])
        for key in (
            "p2698_next_step",
            "fractal_information_ontology",
            "z12_aut_invariant_boundary",
            "premise_based_direction_support",
            "qw2191_still_open",
            "forbidden_promotions",
        ):
            self.assertIn(key, audit["patterns"])
            self.assertGreater(audit["patterns"][key]["count"], 0)

    def test_state_reads_execute_p2698_recommendation_without_replay_promotion(self) -> None:
        state = self.payload["state_reads"]
        self.assertTrue(state["p2698_recommended_p2699"])
        self.assertTrue(state["ontology_candidate_is_internal_only"])
        self.assertTrue(state["existing_direction_support_real_but_premise_based"])
        self.assertTrue(state["h39_qw2191_still_open"])
        self.assertTrue(state["p739_pair12_upgrade_unexported"])
        self.assertTrue(state["p740_pair12_upgrade_unexported"])
        self.assertTrue(state["kappa_is_premise_based_not_used_as_nonpremise_source"])

    def test_aut_z12_calculation_is_finite_and_blocks_directed_generator(self) -> None:
        calc = self.payload["z12_aut_calculation"]
        self.assertEqual(calc["units"], [1, 5, 7, 11])
        self.assertEqual(calc["fixed_points"], [0, 6])
        self.assertEqual(calc["directed_generator_orbit"], [1, 5, 7, 11])
        self.assertTrue(calc["plus_one_and_minus_one_same_orbit"])
        self.assertFalse(calc["can_select_unique_directed_generator_aut_invariantly"])
        self.assertFalse(calc["can_distinguish_plus_one_from_minus_one_aut_invariantly"])
        self.assertEqual(calc["invariant_scalar_basis_dimension"], 6)

    def test_no_go_matrix_and_decision_preserve_no_false_pass(self) -> None:
        rows = self.payload["no_go_matrix"]
        self.assertEqual(len(rows), 4)
        self.assertTrue(all(not row["passes"] for row in rows))
        decision = self.payload["decision"]
        self.assertTrue(decision["bounded_no_go_now"])
        self.assertIn("no-new-live-frontier", decision["next_honest_step"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2699/S1649", MD.read_text(encoding="utf-8"))
        self.assertIn("P2699/S1649", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2699/S1649", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2699/S1649", (ROOT.parent / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
