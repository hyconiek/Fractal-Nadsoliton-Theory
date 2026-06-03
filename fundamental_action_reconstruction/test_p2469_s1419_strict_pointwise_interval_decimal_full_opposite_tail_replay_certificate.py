#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
import unittest
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
import p2469_s1419_strict_pointwise_interval_decimal_full_opposite_tail_replay_certificate as p2469

OUT = ROOT / "generated" / "p2469_s1419_strict_pointwise_interval_decimal_full_opposite_tail_replay_certificate.json"
MD = ROOT / "generated" / "p2469_s1419_strict_pointwise_interval_decimal_full_opposite_tail_replay_certificate.md"


class P2469StrictPointwiseIntervalDecimalFullOppositeTailReplayCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not OUT.exists():
            p2469.main()
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.theorem = cls.payload["strict_pointwise_interval_decimal_full_opposite_tail_replay_certificate"]["theorem_export"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2469")
        self.assertEqual(self.payload["stage_id"], "S1419")
        self.assertIn("FULL_OPPOSITE_TAIL_REPLAY", self.payload["status"])
        self.assertIn("TAIL_ONLY", self.payload["status"])

    def test_inherited_counts(self) -> None:
        self.assertEqual(self.theorem["p2451_total_interval_complement_cell_count_inherited"], 99882)
        self.assertEqual(self.theorem["p2459_total_unreplayed_by_decimal_boundary_chain_cell_count_inherited"], 99846)
        self.assertEqual(self.theorem["p2466_total_tail_boundary_replay_count_inherited"], 6361)
        self.assertEqual(self.theorem["p2467_total_opposite_tail_available_cell_count_inherited"], 45165)
        self.assertEqual(self.theorem["p2467_total_opposite_tail_sentinel_replay_count_inherited"], 42)
        self.assertEqual(self.theorem["p2468_total_opposite_tail_chunked_replay_count_inherited"], 140)
        self.assertEqual(self.theorem["p2468_total_opposite_tail_unreplayed_budget_after_chunked_replay_inherited"], 45025)

    def test_full_opposite_tail_accounting(self) -> None:
        self.assertEqual(self.theorem["total_opposite_tail_full_replay_count"], 45165)
        self.assertEqual(self.theorem["total_opposite_tail_available_cell_count"], 45165)
        self.assertEqual(self.theorem["total_opposite_tail_unreplayed_budget_after_full_tail_replay"], 0)
        self.assertEqual(self.theorem["p2459_unreplayed_budget_not_covered_by_p2466_p2469_tail_replays"], 48320)
        self.assertEqual(self.theorem["p2459_tail_coverage_replayed_cell_count_from_p2466_and_p2469"], 51526)
        self.assertEqual(self.theorem["total_opposite_tail_chunk_count"], 46)

    def test_zero_exclusion_disjointness_and_minimum(self) -> None:
        self.assertTrue(self.theorem["all_opposite_tail_full_replayed_cells_exclude_zero"])
        self.assertGreater(Decimal(self.theorem["minimum_opposite_tail_full_replay_decimal_separation"]), Decimal("0"))
        self.assertTrue(self.theorem["all_opposite_tail_full_replay_indexes_unique_by_family"])
        self.assertTrue(self.theorem["all_opposite_tail_full_replays_disjoint_from_p2466_descent_tail"])
        self.assertTrue(self.theorem["all_family_available_counts_match_p2467"])
        for family in self.theorem["family_full_opposite_tail_replays"]:
            self.assertEqual(family["opposite_tail_full_replay_count"], family["opposite_tail_available_count"])
            self.assertEqual(family["opposite_tail_unreplayed_budget_after_full_tail_replay"], 0)
            self.assertTrue(family["all_opposite_tail_full_replayed_cells_exclude_zero"])
            self.assertGreater(Decimal(family["minimum_opposite_tail_full_replay_decimal_separation"]), Decimal("0"))
            self.assertTrue(family["full_opposite_tail_replay_exported_by_this_family"])

    def test_gatekeepers_and_negative_controls(self) -> None:
        for name, value in self.payload["gatekeeper_checks"].items():
            self.assertTrue(value, name)
        self.assertTrue(self.theorem["full_opposite_tail_replay_exported_by_this_certificate"])
        self.assertFalse(self.theorem["remaining_complement_segments_replayed_by_this_certificate"])
        self.assertFalse(self.theorem["decimal_full_complement_replay_exported_by_this_certificate"])
        self.assertFalse(self.theorem["directed_rounding_interval_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["symbolic_root_exclusion_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["pointwise_coordinate_selector_exported_by_this_certificate"])
        self.assertFalse(self.theorem["strict_observable_source_constraint_exported_by_this_certificate"])
        self.assertFalse(self.theorem["gauge_slice_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["strict_physical_value_generator_exported"])
        self.assertFalse(self.theorem["qw2191_discharged"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["legacy_role_transfer_exported"])
        self.assertFalse(self.theorem["toe_closure_exported"])

    def test_rg_audit_and_docs(self) -> None:
        audit = self.payload["rg_non_duplication_audit"]
        self.assertEqual(audit["tool"], "rg")
        self.assertIn("new_packet", audit["patterns"])
        self.assertIn("full opposite-tail replay", MD.read_text(encoding="utf-8"))
        self.assertIn("P2469/S1419", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2469/S1419", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
