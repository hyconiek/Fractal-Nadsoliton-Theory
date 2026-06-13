from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2694_s1644_fresh_broad_state_map_selection_after_p2680_closure.py"
OUT = ROOT / "generated" / "p2694_s1644_fresh_broad_state_map_selection_after_p2680_closure.json"
MD = ROOT / "generated" / "p2694_s1644_fresh_broad_state_map_selection_after_p2680_closure.md"


class P2694FreshBroadStateMapSelectionAfterP2680ClosureTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_grep_covers_broad_state_map(self) -> None:
        audit = self.payload["content_grep"]
        self.assertIn("fresh broad state-map", audit["mode"])
        for key in (
            "p2693_selected_p2694",
            "f3_residual_direct_g_family",
            "closed_bridge_source_round",
            "closed_lagrangian_lower_selector",
            "forbidden_reopen",
        ):
            self.assertIn(key, audit["patterns"])
            self.assertGreater(audit["patterns"][key]["count"], 0)

    def test_state_reads_support_fresh_selection(self) -> None:
        state = self.payload["state_reads"]
        self.assertTrue(state["p2693_selected_p2694"])
        self.assertTrue(state["p2693_only_fresh_state_map_live"])
        self.assertTrue(state["p2682_m2_target_bounded_no_go"])
        self.assertFalse(state["p2684_lower_boundary_primary"])
        self.assertTrue(state["p2687_freezes_lagrangian_lane"])
        self.assertTrue(state["p2688_selected_bridge_uv_lane"])
        self.assertTrue(state["f3_names_g4_g6_gY_zero_witnesses"])

    def test_residual_direct_target_matrix_selects_three_non_m2_targets(self) -> None:
        rows = self.payload["residual_direct_target_matrix"]
        self.assertEqual(len(rows), 3)
        self.assertTrue(all(row["live_now"] for row in rows))
        self.assertTrue(all(row["finite_computable"] for row in rows))
        self.assertTrue(all(not row["already_attacked_by_ax13_p631_m2_freeze"] for row in rows))
        self.assertEqual(
            {row["target_id"] for row in rows},
            {
                "direct_g4_c1s1_shift_defect_zero_witness",
                "direct_g6_c1s1_shift_defect_zero_witness",
                "direct_gY_c1s1_shift_defect_zero_witness",
            },
        )

    def test_broad_lane_matrix_keeps_only_residual_direct_lane_live(self) -> None:
        lanes = {row["lane"]: row for row in self.payload["broad_lane_matrix"]}
        self.assertFalse(lanes["generic_bridge_or_p2680_nonselector_replay"]["live_now"])
        self.assertFalse(lanes["direct_m2_psi4_attacked_target"]["live_now"])
        self.assertTrue(lanes["direct_residual_g_family_zero_witnesses"]["live_now"])
        self.assertTrue(lanes["direct_residual_g_family_zero_witnesses"]["proof_grade_next"])
        self.assertFalse(lanes["lagrangian_eom_reverse_closure"]["live_now"])
        self.assertFalse(lanes["lower_boundary_pair12_cycle"]["live_now"])
        self.assertFalse(lanes["selector_tau_pair_or_beta_tors_replay"]["live_now"])
        self.assertFalse(lanes["role_transfer_ltotal_toe"]["live_now"])

    def test_selection_recommends_p2695_without_forbidden_reopens(self) -> None:
        selection = self.payload["selection"]
        self.assertEqual(selection["selected_lane"], "direct_residual_g_family_zero_witnesses")
        self.assertEqual(selection["selected_next_packet"], "P2695_residual_direct_g4_g6_gY_c1s1_zero_witness_no_go_matrix")
        self.assertEqual(len(selection["selected_targets"]), 3)
        self.assertIn("P2695", selection["next_honest_step"])
        self.assertFalse(selection["bridge_reopened_now"])
        self.assertFalse(selection["m2_psi4_reopened_now"])
        self.assertFalse(selection["selector_replay_now"])
        self.assertFalse(selection["ltotal_promoted_now"])
        self.assertFalse(selection["toe_closed_now"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2694/S1644", MD.read_text(encoding="utf-8"))
        self.assertIn("P2694/S1644", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2694/S1644", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2694/S1644", (ROOT.parent / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
