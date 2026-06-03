from __future__ import annotations

import json
import sys
import unittest
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2481_s1431_strict_pointwise_interval_decimal_p2459_boundary_handoff_collar_replay_audit import (
    MD,
    OUT,
    build_payload,
    write_markdown,
)


class P2481StrictPointwiseIntervalDecimalP2459BoundaryHandoffCollarReplayAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["strict_pointwise_interval_decimal_p2459_boundary_handoff_collar_replay_audit"]["theorem_export"]
        cls.replay = cls.theorem["handoff_collar_replay"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2481")
        self.assertEqual(self.payload["stage_id"], "S1431")
        self.assertIn("NO_FULL_COMPLEMENT", self.payload["status"])
        self.assertIn("NO_DIRECTED_ROUNDING", self.payload["status"])

    def test_handoff_collar_counts_and_binding(self) -> None:
        self.assertEqual(self.theorem["parent_segment_index"], 2)
        self.assertEqual(self.theorem["parent_uncovered_index"], 0)
        self.assertEqual(self.theorem["fresh_boundary_band_cell_replay_count"], 6)
        self.assertEqual(self.theorem["fresh_dyadic_subcell_replay_count"], 128)
        self.assertEqual(self.theorem["total_fresh_handoff_collar_replay_rows"], 134)
        self.assertEqual(Decimal(self.theorem["handoff_endpoint_gap_boundary_to_parent"]), Decimal("0"))
        self.assertTrue(self.theorem["p2456_boundary_band_cells_are_inherited_covered_boundary_chain_not_new_p2459_unreplayed_cells"])

    def test_zero_exclusion_and_monotone_collar_order(self) -> None:
        self.assertTrue(self.theorem["all_collar_rows_exclude_zero"])
        self.assertTrue(self.theorem["all_collar_rows_have_positive_separation"])
        self.assertTrue(self.theorem["all_consecutive_collar_separations_strictly_increase_left_to_right"])
        self.assertTrue(self.theorem["all_consecutive_collar_rows_are_exactly_adjacent"])
        self.assertEqual(Decimal(self.theorem["minimum_consecutive_endpoint_gap"]), Decimal("0"))
        self.assertEqual(Decimal(self.theorem["maximum_consecutive_endpoint_gap"]), Decimal("0"))
        self.assertGreater(Decimal(self.theorem["minimum_consecutive_positive_delta"]), Decimal("0"))
        self.assertGreater(Decimal(self.theorem["minimum_collar_decimal_separation"]), Decimal("0"))
        self.assertTrue(self.theorem["minimum_is_p2456_right_boundary_band_leftmost_cell"])
        self.assertTrue(self.theorem["maximum_is_p2480_rightmost_dyadic_subcell"])
        self.assertEqual(self.replay["minimum_collar_replay_row"]["collar_side"], "p2456_right_boundary_band")

    def test_coverage_budget_and_inherited_counts(self) -> None:
        self.assertEqual(self.theorem["p2471_p2459_universe_count_inherited"], 99846)
        self.assertEqual(self.theorem["p2481_fresh_decimal_evaluation_row_count_not_a_coverage_count"], 134)
        self.assertEqual(self.theorem["p2481_fresh_decimal_evaluation_row_ratio_not_a_p2459_coverage_fraction"], "134/99846")
        self.assertEqual(self.theorem["targeted_p2481_new_p2459_unreplayed_parent_cell_scope_against_p2459_universe"], "1/99846")
        self.assertEqual(self.theorem["targeted_p2481_new_p2459_unreplayed_parent_cell_count"], 1)
        self.assertEqual(self.theorem["p2481_subcell_rows_inside_single_parent_cell_not_distinct_p2459_cells"], 128)
        self.assertEqual(self.theorem["p2456_inherited_covered_boundary_chain_cell_count"], 6)
        self.assertEqual(self.theorem["p2480_dyadic_subcell_count_inherited"], 128)
        self.assertEqual(self.theorem["p2479_prefix_replay_count_inherited"], 1116)
        self.assertEqual(self.theorem["p2479_prefix_plus_p2478_union_count_inherited"], 1132)
        self.assertTrue(self.theorem["no_full_complement_claimed_by_this_certificate"])
        self.assertFalse(self.theorem["full_complement_replay_exported_by_this_certificate"])

    def test_gatekeepers_and_negative_controls(self) -> None:
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertFalse(self.theorem["directed_rounding_interval_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["symbolic_root_exclusion_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["analytic_monotonicity_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["global_continuum_root_exclusion_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["pointwise_coordinate_selector_exported_by_this_certificate"])
        self.assertFalse(self.theorem["strict_observable_source_constraint_exported_by_this_certificate"])
        self.assertFalse(self.theorem["gauge_slice_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["strict_physical_value_generator_exported"])
        self.assertFalse(self.theorem["qw2191_discharged"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["legacy_role_transfer_exported"])
        self.assertFalse(self.theorem["toe_closure_exported"])

    def test_rg_audit_and_documentation_markers(self) -> None:
        audit = self.payload["rg_non_duplication_audit"]
        self.assertEqual(audit["tool"], "rg")
        self.assertIn("new_packet", audit["patterns"])
        self.assertGreaterEqual(audit["patterns"]["new_packet"]["count"], 0)
        self.assertIn("P2481/S1431", MD.read_text(encoding="utf-8"))
        self.assertIn("P2481/S1431", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2481/S1431", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
