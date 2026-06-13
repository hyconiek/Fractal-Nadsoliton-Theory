from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2697_s1647_post_direct_route_state_map_no_new_live_frontier_reconciliation.py"
OUT = ROOT / "generated" / "p2697_s1647_post_direct_route_state_map_no_new_live_frontier_reconciliation.json"
MD = ROOT / "generated" / "p2697_s1647_post_direct_route_state_map_no_new_live_frontier_reconciliation.md"


class P2697PostDirectRouteStateMapNoNewLiveFrontierTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_grep_covers_broad_state_map(self) -> None:
        audit = self.payload["content_grep"]
        self.assertIn("no-new-live-frontier", audit["mode"])
        for key in ("p2696_selected_p2697", "direct_route_closed", "closed_nonselector_bridge", "closed_lagrangian_lower_boundary", "selector_h37_not_new", "forbidden_promotions"):
            self.assertIn(key, audit["patterns"])
            self.assertGreater(audit["patterns"][key]["count"], 0)

    def test_state_reads_support_no_new_frontier(self) -> None:
        state = self.payload["state_reads"]
        self.assertTrue(state["p2696_selected_p2697"])
        self.assertTrue(state["direct_g_family_closed_by_p2695"])
        self.assertTrue(state["pair1_c1c1_s1s1_bounded_no_go_by_p2696"])
        self.assertTrue(state["p631_direct_formal_negative_freeze"])
        self.assertTrue(state["p2693_nonselector_source_inventory_closed"])
        self.assertTrue(state["p2687_lagrangian_frozen"])
        self.assertTrue(state["p2684_lower_boundary_not_primary"])
        self.assertTrue(state["p2683_h37_not_missing_replay"])
        self.assertFalse(state["p739_strict_core_upgrade_exported"])
        self.assertFalse(state["p740_strict_core_upgrade_exported"])

    def test_freshness_gate_matrix_has_no_live_lane(self) -> None:
        rows = self.payload["freshness_gate_matrix"]
        self.assertGreaterEqual(len(rows), 6)
        self.assertTrue(all(row["closed_or_replay_gated_now"] for row in rows))
        self.assertTrue(all(not row["live_now"] for row in rows))
        lanes = {row["lane"] for row in rows}
        self.assertIn("direct_route_g_family_and_pair1_residual_cancellation", lanes)
        self.assertIn("selector_h37_t171_qw2191_replay", lanes)
        self.assertIn("role_transfer_ltotal_toe", lanes)

    def test_decision_exports_no_new_live_frontier_certificate(self) -> None:
        decision = self.payload["decision"]
        self.assertTrue(decision["no_new_live_frontier_certificate"])
        self.assertEqual(decision["live_lanes_now"], [])
        self.assertIn("genuinely new strict typed provider", decision["next_honest_step"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2697/S1647", MD.read_text(encoding="utf-8"))
        self.assertIn("P2697/S1647", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2697/S1647", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2697/S1647", (ROOT.parent / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
