from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2620_s1570_p2618_p2619_bridge_two_obstruction_cut import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2620BridgeTwoObstructionCutTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["p2618_p2619_bridge_two_obstruction_cut"]["theorem_export"]

    def test_identity_and_inherited_blocks(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2620")
        self.assertEqual(self.payload["stage_id"], "S1570")
        self.assertIn("TWO_OBSTRUCTION_BRIDGE_SOURCE_CUT", self.payload["status"])
        self.assertTrue(self.theorem["inherits_p2618_full_bridge_block"])
        self.assertTrue(self.theorem["inherits_p2619_selector_source_block"])
        self.assertTrue(self.theorem["inherits_p2616_role_obstruction"])

    def test_truth_table_has_exactly_one_accepting_row(self) -> None:
        self.assertEqual(self.theorem["bridge_cut_truth_table_row_count"], 4)
        self.assertEqual(self.theorem["bridge_cut_accepting_row_count"], 1)
        accepting = self.theorem["bridge_cut_accepting_row"]
        self.assertTrue(accepting["assignment"]["nonlinear_damping_completion_source"])
        self.assertTrue(accepting["assignment"]["orientation_odd_selector_source"])
        self.assertFalse(self.theorem["current_bridge_source_cut_repaired"])

    def test_shortcuts_rejected_except_both_sources(self) -> None:
        rows = {row["shortcut"]: row for row in self.theorem["shortcut_rejection_rows"]}
        self.assertFalse(rows["eta_9_5_exponent_source_only"]["passes_bridge_source_cut"])
        self.assertFalse(rows["beta_tors_scalar_renormalization"]["passes_bridge_source_cut"])
        self.assertFalse(rows["axis_only_selector_up_to_Z2"]["passes_bridge_source_cut"])
        self.assertFalse(rows["GF2_bookkeeping_rank_or_cycle_data"]["passes_bridge_source_cut"])
        self.assertFalse(rows["damping_source_plus_no_selector"]["passes_bridge_source_cut"])
        self.assertFalse(rows["selector_source_plus_no_damping_completion"]["passes_bridge_source_cut"])
        self.assertTrue(rows["both_independent_sources_supplied"]["passes_bridge_source_cut"])

    def test_scope_guards(self) -> None:
        self.assertFalse(self.theorem["full_bridge_revalidated"])
        self.assertFalse(self.theorem["role_transfer_revalidated"])
        self.assertFalse(self.theorem["role_bearing_ltotal_reenabled"])
        self.assertFalse(self.theorem["beta_tors_chi11_route_reopened"])
        self.assertFalse(self.theorem["gf2_bridge_revalidated"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_packet"])
        self.assertFalse(self.theorem["toe_closure_claimed"])

    def test_gatekeepers_and_docs(self) -> None:
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2620/S1570", MD.read_text(encoding="utf-8"))
        self.assertIn("P2620/S1570", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2620/S1570", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
