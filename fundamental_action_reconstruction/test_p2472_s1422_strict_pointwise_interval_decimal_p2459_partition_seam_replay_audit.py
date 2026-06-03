#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
import unittest
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
import p2472_s1422_strict_pointwise_interval_decimal_p2459_partition_seam_replay_audit as p2472

OUT = ROOT / "generated" / "p2472_s1422_strict_pointwise_interval_decimal_p2459_partition_seam_replay_audit.json"
MD = ROOT / "generated" / "p2472_s1422_strict_pointwise_interval_decimal_p2459_partition_seam_replay_audit.md"


class P2472StrictPointwiseIntervalDecimalP2459PartitionSeamReplayAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not OUT.exists():
            p2472.main()
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.theorem = cls.payload["strict_pointwise_interval_decimal_p2459_partition_seam_replay_audit"]["theorem_export"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2472")
        self.assertEqual(self.payload["stage_id"], "S1422")
        self.assertIn("P2459_PARTITION_SEAM_REPLAY_AUDIT", self.payload["status"])
        self.assertIn("NO_DIRECTED_ROUNDING", self.payload["status"])

    def test_inherited_counts_and_seam_totals(self) -> None:
        self.assertEqual(self.theorem["p2459_total_unreplayed_by_decimal_boundary_chain_cell_count_inherited"], 99846)
        self.assertEqual(self.theorem["p2466_total_tail_boundary_replay_count_inherited"], 6361)
        self.assertEqual(self.theorem["p2469_total_opposite_tail_full_replay_count_inherited"], 45165)
        self.assertEqual(self.theorem["p2470_total_remaining_non_tail_full_replay_count_inherited"], 48320)
        self.assertEqual(self.theorem["p2471_total_p2459_universe_cell_count_rebuilt_inherited"], 99846)
        self.assertTrue(self.theorem["p2471_all_family_partitions_are_disjoint_and_complete_inherited"])
        self.assertEqual(self.theorem["total_transition_pair_count"], 4)
        self.assertEqual(self.theorem["total_partition_seam_replay_count"], 16)

    def test_transition_adjacency_and_zero_exclusion(self) -> None:
        self.assertTrue(self.theorem["all_transition_pairs_are_adjacent"])
        self.assertTrue(self.theorem["all_seam_replayed_cells_exclude_zero"])
        self.assertGreater(Decimal(self.theorem["minimum_partition_seam_decimal_separation"]), Decimal("0"))
        expected_counts = {
            "zero_projection_amplitude": (2, 9),
            "stationary_factor": (2, 7),
        }
        for family in self.theorem["family_partition_seam_replays"]:
            self.assertEqual((family["transition_pair_count"], family["seam_replay_count"]), expected_counts[family["family"]])
            self.assertTrue(family["all_transition_pairs_are_adjacent"])
            self.assertTrue(family["all_seam_replayed_cells_exclude_zero"])
            self.assertGreater(Decimal(family["minimum_seam_decimal_separation"]), Decimal("0"))
            for pair in family["transition_pairs"]:
                self.assertTrue(pair["adjacent_indexes"])

    def test_gatekeepers_and_negative_controls(self) -> None:
        for name, value in self.payload["gatekeeper_checks"].items():
            self.assertTrue(value, name)
        self.assertTrue(self.theorem["finite_partition_seam_replay_audit_exported"])
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

    def test_rg_audit_docs_and_lay_summary(self) -> None:
        audit = self.payload["rg_non_duplication_audit"]
        self.assertEqual(audit["tool"], "rg")
        self.assertIn("new_packet", audit["patterns"])
        md_text = MD.read_text(encoding="utf-8")
        self.assertIn("Plain-language progress note", md_text)
        self.assertIn("off-by-one seam error", self.theorem["lay_summary"])
        self.assertIn("P2472/S1422", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2472/S1422", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
