from __future__ import annotations

import json
import sys
import unittest
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2484_s1434_strict_pointwise_interval_decimal_p2459_boundary_side_dyadic_secant_margin_replay_audit import (
    MD,
    OUT,
    build_payload,
    write_markdown,
)


class P2484StrictPointwiseIntervalDecimalP2459BoundarySideDyadicSecantMarginReplayAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["strict_pointwise_interval_decimal_p2459_boundary_side_dyadic_secant_margin_replay_audit"]["theorem_export"]
        cls.replay = cls.theorem["boundary_side_dyadic_secant_margin_replay"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2484")
        self.assertEqual(self.payload["stage_id"], "S1434")
        self.assertIn("NO_ROOT_WINDOW", self.payload["status"])
        self.assertIn("NO_ANALYTIC_MONOTONICITY", self.payload["status"])
        self.assertIn("NO_COVERAGE_INCREASE", self.payload["status"])
        self.assertIn("NO_DIRECTED_ROUNDING", self.payload["status"])

    def test_parent_binding_and_secant_levels(self) -> None:
        self.assertEqual(self.theorem["parent_collar_side"], "p2456_right_boundary_band")
        self.assertEqual(self.theorem["parent_boundary_band_index"], 0)
        self.assertEqual(self.theorem["parent_dyadic_subcell_index"], 0)
        self.assertEqual(self.theorem["secant_level_count"], 32)
        self.assertEqual(self.theorem["consecutive_secant_pair_count"], 31)
        self.assertEqual(self.theorem["tightest_secant_level"], 32)
        self.assertEqual(self.theorem["tightest_width_fraction_of_p2482_leftmost_subcell"], "1/4294967296")
        self.assertTrue(self.theorem["weakest_secant_row_is_coarsest_level"])
        self.assertEqual(self.theorem["p2483_nested_level_count_inherited"], 16)

    def test_zero_exclusion_width_halving_and_secant_margins(self) -> None:
        self.assertTrue(self.theorem["all_secant_rows_exclude_zero"])
        self.assertTrue(self.theorem["all_secant_rows_have_positive_separation"])
        self.assertTrue(self.theorem["all_secant_rows_share_left_boundary_anchor"])
        self.assertTrue(self.theorem["all_consecutive_widths_halve"])
        self.assertTrue(self.theorem["all_consecutive_lower_bounds_strictly_increase"])
        self.assertTrue(self.theorem["all_consecutive_secant_margins_positive"])
        self.assertGreater(Decimal(self.theorem["minimum_consecutive_positive_lower_bound_delta"]), Decimal("0"))
        self.assertGreater(Decimal(self.theorem["minimum_positive_secant_margin_per_removed_width"]), Decimal("0"))
        self.assertGreater(Decimal(self.theorem["maximum_positive_secant_margin_per_removed_width"]), Decimal(self.theorem["minimum_positive_secant_margin_per_removed_width"]))
        self.assertGreater(Decimal(self.theorem["weakest_secant_decimal_separation"]), Decimal("0"))
        self.assertGreater(Decimal(self.theorem["tightest_secant_decimal_separation"]), Decimal(self.theorem["weakest_secant_decimal_separation"]))
        self.assertTrue(self.theorem["tightest_lower_bound_exceeds_p2482_leftmost_subcell_bound"])
        self.assertGreater(Decimal(self.theorem["tightest_minus_p2482_leftmost_subcell_lower_bound_delta"]), Decimal("0"))

    def test_coverage_budget_and_inherited_counts(self) -> None:
        self.assertEqual(self.theorem["p2471_p2459_universe_count_inherited"], 99846)
        self.assertEqual(self.theorem["p2484_fresh_decimal_evaluation_row_count_not_a_coverage_count"], 32)
        self.assertEqual(self.theorem["p2484_consecutive_secant_margin_count_not_a_coverage_count"], 31)
        self.assertEqual(self.theorem["p2484_fresh_decimal_evaluation_row_ratio_not_a_p2459_coverage_fraction"], "32/99846")
        self.assertEqual(self.theorem["targeted_p2484_new_p2459_unreplayed_cell_count"], 0)
        self.assertEqual(self.theorem["targeted_p2484_new_p2459_unreplayed_cell_scope_against_p2459_universe"], "0/99846")
        self.assertTrue(self.theorem["p2484_refines_one_inherited_p2456_covered_boundary_chain_cell"])
        self.assertTrue(self.theorem["no_full_complement_claimed_by_this_certificate"])
        self.assertFalse(self.theorem["full_complement_replay_exported_by_this_certificate"])

    def test_gatekeepers_and_negative_controls(self) -> None:
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertFalse(self.theorem["directed_rounding_interval_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["symbolic_root_exclusion_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["analytic_monotonicity_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["root_window_exclusion_theorem_exported_by_this_certificate"])
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
        self.assertIn("P2484/S1434", MD.read_text(encoding="utf-8"))
        self.assertIn("P2484/S1434", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2484/S1434", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
