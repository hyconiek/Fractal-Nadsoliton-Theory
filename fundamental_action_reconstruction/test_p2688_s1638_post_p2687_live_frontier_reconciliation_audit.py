from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2688_s1638_post_p2687_live_frontier_reconciliation_audit.py"
OUT = ROOT / "generated" / "p2688_s1638_post_p2687_live_frontier_reconciliation_audit.json"
MD = ROOT / "generated" / "p2688_s1638_post_p2687_live_frontier_reconciliation_audit.md"


class P2688PostP2687LiveFrontierReconciliationAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_grep_covers_reconciled_state_map(self) -> None:
        audit = self.payload["content_grep"]
        self.assertIn("post-P2687 live-frontier reconciliation", audit["mode"])
        for key in (
            "post_p2687_state_map",
            "direct_route_status",
            "lagrangian_eom_freeze",
            "lower_boundary_freeze",
            "non_selector_bridge_atoms",
            "forbidden_replays",
        ):
            self.assertIn(key, audit["patterns"])
            self.assertGreater(audit["patterns"][key]["count"], 0)

    def test_state_reads_override_stale_p46_return(self) -> None:
        state = self.payload["state_reads"]
        self.assertTrue(state["p2687_freezes_lagrangian_lane"])
        self.assertTrue(state["p2687_recommended_stale_p46_pivot"])
        self.assertTrue(state["ax13_target_eom_closed_external"])
        self.assertFalse(state["p51_direct_route_full_closure_pass"])
        self.assertTrue(state["p631_direct_route_negative_freeze"])
        self.assertTrue(state["p51_remaining_g_family_witnesses_absent"])

    def test_bridge_uv_unit_atom_is_selected_live_frontier(self) -> None:
        state = self.payload["state_reads"]
        self.assertIn("canonical_length_or_uv_unit_source", state["p2680_missing_non_selector_atoms"])
        self.assertFalse(state["p2650_exports_canonical_unit"])
        self.assertFalse(state["p2653_exports_typed_metric_uv_source"])
        self.assertTrue(state["p2662_conditional_unit_map_present"])
        self.assertTrue(state["p2663_one_bit_carrier_present"])
        rows = {row["lane_id"]: row for row in self.payload["lane_matrix"]}
        self.assertTrue(rows["non_selector_bridge_canonical_uv_unit_atom"]["live_now"])
        self.assertTrue(rows["non_selector_bridge_canonical_uv_unit_atom"]["proof_grade_next"])
        self.assertFalse(rows["direct_p46_p50_m2_psi4"]["live_now"])
        self.assertFalse(rows["lagrangian_eom_anisotropic_reverse_closure"]["live_now"])

    def test_selection_is_p2689_without_closure_exports(self) -> None:
        selection = self.payload["selection"]
        self.assertEqual(selection["selected_lane"], "non_selector_bridge_canonical_uv_unit_atom")
        self.assertIn("P2689", selection["selected_next_packet"])
        self.assertTrue(selection["stale_p46_recommendation_overridden"])
        self.assertFalse(selection["uv_unit_or_beta_source_exported_now"])
        self.assertFalse(selection["role_transfer_started_now"])
        self.assertFalse(selection["toe_closed_now"])
        self.assertTrue(all(not item["satisfied_now"] for item in selection["selected_obligations"]))
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))

    def test_documents_updated(self) -> None:
        self.assertIn("P2688/S1638", MD.read_text(encoding="utf-8"))
        self.assertIn("P2688/S1638", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2688/S1638", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2688/S1638", (ROOT.parent / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
