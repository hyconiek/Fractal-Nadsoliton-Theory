from __future__ import annotations

import json
import sys
import unittest
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2477_s1427_strict_pointwise_interval_decimal_p2459_lower_neighbor_exception_expanded_halo_replay_audit import (
    MD,
    OUT,
    build_payload,
    write_markdown,
)


class P2477StrictPointwiseIntervalDecimalP2459LowerNeighborExceptionExpandedHaloReplayAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["strict_pointwise_interval_decimal_p2459_lower_neighbor_exception_expanded_halo_replay_audit"]["theorem_export"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2477")
        self.assertEqual(self.payload["stage_id"], "S1427")
        self.assertIn("NO_FULL_COMPLEMENT", self.payload["status"])
        self.assertIn("NO_DIRECTED_ROUNDING", self.payload["status"])

    def test_expanded_exception_replay_counts_and_minimum(self) -> None:
        self.assertEqual(self.theorem["lower_neighbor_exception_count_inherited_from_p2476"], 1)
        self.assertEqual(self.theorem["expanded_radius"], 4)
        self.assertEqual(self.theorem["total_unique_expanded_exception_replay_count"], 11)
        self.assertEqual(self.theorem["incremental_cells_beyond_p2475_radius_two_halo_count"], 6)
        self.assertGreater(Decimal(self.theorem["minimum_expanded_exception_decimal_separation"]), Decimal("0"))
        replay = self.theorem["expanded_exception_replays"][0]
        self.assertEqual(replay["minimum_expanded_exception_replay_row"]["uncovered_index"], 1131)
        self.assertTrue(replay["lowest_is_left_boundary_of_expanded_replay"])
        self.assertTrue(replay["center_has_lower_left_flank_in_expanded_replay"])
        self.assertFalse(replay["center_has_lower_right_flank_in_expanded_replay"])
        self.assertTrue(replay["all_consecutive_pairs_strictly_increase_left_to_right"])

    def test_zero_exclusion_and_coverage_budget(self) -> None:
        self.assertTrue(self.theorem["all_expanded_exception_cells_exclude_zero"])
        self.assertTrue(self.theorem["all_expanded_exception_cells_have_positive_separation"])
        self.assertEqual(self.theorem["p2471_p2459_universe_count_inherited"], 99846)
        self.assertEqual(self.theorem["targeted_p2477_replay_fraction_of_p2459_universe"], "11/99846")
        self.assertEqual(self.theorem["targeted_p2477_residual_not_freshly_replayed_count_against_p2459_universe"], 99835)
        self.assertEqual(self.theorem["finite_chain_coverage_budget_inherited_from_p2471"], 0)
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
        self.assertIn("P2477/S1427", MD.read_text(encoding="utf-8"))
        self.assertIn("P2477/S1427", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2477/S1427", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
