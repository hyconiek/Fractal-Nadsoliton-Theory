from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2693_s1643_post_p2680_nonselector_source_inventory_closure_state_map.py"
OUT = ROOT / "generated" / "p2693_s1643_post_p2680_nonselector_source_inventory_closure_state_map.json"
MD = ROOT / "generated" / "p2693_s1643_post_p2680_nonselector_source_inventory_closure_state_map.md"


class P2693PostP2680NonselectorSourceInventoryClosureStateMapTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_grep_covers_inventory_closure(self) -> None:
        audit = self.payload["content_grep"]
        self.assertIn("post-P2680", audit["mode"])
        for key in (
            "p2692_selected_p2693",
            "p2680_nonselector_atoms",
            "uv_unit_route_freeze",
            "amplitude_beta_freeze",
            "forbidden_replay_imports",
        ):
            self.assertIn(key, audit["patterns"])
            self.assertGreater(audit["patterns"][key]["count"], 0)

    def test_state_reads_show_all_named_nonselector_atoms_unexported(self) -> None:
        state = self.payload["state_reads"]
        self.assertTrue(state["p2692_selected_p2693"])
        self.assertFalse(state["p2680_amplitude_source_exported"])
        self.assertFalse(state["p2680_canonical_uv_source_exported"])
        self.assertFalse(state["p2680_positive_beta_source_exported"])
        self.assertFalse(state["p2689_uv_unit_or_beta_source_exported"])
        self.assertTrue(state["p2690_freezes_entropy_uv_route"])
        self.assertTrue(state["p2691_bounds_alpha_atom"])
        self.assertTrue(state["p2692_bounds_beta_atom"])
        self.assertFalse(state["p2692_beta_source_exported"])

    def test_nonselector_atom_status_marks_three_bounded_no_go_rows(self) -> None:
        rows = self.payload["nonselector_atom_status"]
        self.assertEqual(len(rows), 3)
        self.assertTrue(all(row["formal_near_miss_present"] for row in rows))
        self.assertTrue(all(row["bounded_no_go_now"] for row in rows))
        self.assertFalse(any(row["source_exported_now"] for row in rows))
        self.assertEqual(
            {row["atom"] for row in rows},
            {
                "canonical_length_or_uv_unit_source",
                "amplitude_role_safe_source",
                "target_independent_positive_beta_or_z_beta_source",
            },
        )

    def test_bridge_source_lattice_requires_three_new_true_bits(self) -> None:
        lattice = self.payload["bridge_source_lattice"]
        self.assertEqual(lattice["total_states"], 8)
        self.assertEqual(lattice["passing_states"], 1)
        self.assertEqual(lattice["hamming_distance_to_nonselector_bridge_source_pass"], 3)
        self.assertFalse(lattice["nonselector_bridge_source_gate_passes_now"])
        self.assertTrue(all(value is False for value in lattice["current_state"].values()))

    def test_lane_reconciliation_allows_only_fresh_broad_state_map(self) -> None:
        lanes = {row["lane"]: row for row in self.payload["lane_reconciliation"]}
        self.assertFalse(lanes["generic_legacy_to_strict_bridge_completion"]["live_now"])
        self.assertFalse(lanes["p2680_nonselector_bridge_source_atoms"]["live_now"])
        self.assertFalse(lanes["selector_or_tau_pair_replay"]["live_now"])
        self.assertFalse(lanes["role_transfer_or_ltotal_promotion"]["live_now"])
        self.assertTrue(lanes["fresh_broad_state_map_selection"]["live_now"])

    def test_decision_recommends_p2694_and_updates_docs(self) -> None:
        decision = self.payload["decision"]
        self.assertFalse(decision["nonselector_bridge_source_gate_passes_now"])
        self.assertEqual(decision["hamming_distance_to_nonselector_bridge_source_pass"], 3)
        self.assertEqual(decision["live_lanes"], ["fresh_broad_state_map_selection"])
        self.assertIn("P2694", decision["next_honest_step"])
        self.assertFalse(decision["bridge_completion_claimed_now"])
        self.assertFalse(decision["role_transfer_started_now"])
        self.assertFalse(decision["ltotal_promoted_now"])
        self.assertFalse(decision["toe_closed_now"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2693/S1643", MD.read_text(encoding="utf-8"))
        self.assertIn("P2693/S1643", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2693/S1643", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2693/S1643", (ROOT.parent / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
