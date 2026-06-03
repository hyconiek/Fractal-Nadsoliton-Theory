from __future__ import annotations

import json
import sys
import unittest
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2480_s1430_strict_pointwise_interval_decimal_p2459_segment_start_cell_dyadic_refinement_replay_audit import (
    MD,
    OUT,
    build_payload,
    write_markdown,
)


class P2480StrictPointwiseIntervalDecimalP2459SegmentStartCellDyadicRefinementReplayAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["strict_pointwise_interval_decimal_p2459_segment_start_cell_dyadic_refinement_replay_audit"]["theorem_export"]
        cls.replay = cls.theorem["dyadic_refinement_replay"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2480")
        self.assertEqual(self.payload["stage_id"], "S1430")
        self.assertIn("NO_FULL_COMPLEMENT", self.payload["status"])
        self.assertIn("NO_DIRECTED_ROUNDING", self.payload["status"])

    def test_dyadic_refinement_counts_and_order(self) -> None:
        self.assertEqual(self.theorem["dyadic_subcell_count"], 128)
        self.assertEqual(self.theorem["parent_segment_index"], 2)
        self.assertEqual(self.theorem["parent_uncovered_index"], 0)
        self.assertEqual(self.replay["minimum_subcell_replay_row"]["dyadic_subcell_index"], 0)
        self.assertEqual(self.replay["maximum_subcell_replay_row"]["dyadic_subcell_index"], 127)
        self.assertTrue(self.theorem["minimum_is_leftmost_subcell"])
        self.assertTrue(self.theorem["maximum_is_rightmost_subcell"])
        self.assertTrue(self.theorem["all_consecutive_subcell_separations_strictly_increase_left_to_right"])
        self.assertGreater(Decimal(self.theorem["minimum_consecutive_positive_delta"]), Decimal("0"))
        self.assertGreater(Decimal(self.theorem["minimum_subcell_decimal_separation"]), Decimal("0"))

    def test_zero_exclusion_and_coverage_budget(self) -> None:
        self.assertTrue(self.theorem["all_subcells_exclude_zero"])
        self.assertTrue(self.theorem["all_subcells_have_positive_separation"])
        self.assertEqual(self.theorem["p2471_p2459_universe_count_inherited"], 99846)
        self.assertEqual(self.theorem["targeted_p2480_subcell_replay_fraction_against_p2459_universe"], "128/99846")
        self.assertEqual(self.theorem["targeted_p2480_parent_cell_count_against_p2459_universe"], "1/99846")
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
        self.assertIn("P2480/S1430", MD.read_text(encoding="utf-8"))
        self.assertIn("P2480/S1430", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2480/S1430", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
