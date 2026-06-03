#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
import unittest
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
import p2475_s1425_strict_pointwise_interval_decimal_p2459_critical_minimum_halo_replay_audit as p2475

OUT = ROOT / "generated" / "p2475_s1425_strict_pointwise_interval_decimal_p2459_critical_minimum_halo_replay_audit.json"
MD = ROOT / "generated" / "p2475_s1425_strict_pointwise_interval_decimal_p2459_critical_minimum_halo_replay_audit.md"


class P2475StrictPointwiseIntervalDecimalP2459CriticalMinimumHaloReplayAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not OUT.exists():
            p2475.main()
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.theorem = cls.payload["strict_pointwise_interval_decimal_p2459_critical_minimum_halo_replay_audit"]["theorem_export"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2475")
        self.assertEqual(self.payload["stage_id"], "S1425")
        self.assertIn("P2459_CRITICAL_MINIMUM_HALO_REPLAY_AUDIT", self.payload["status"])
        self.assertIn("NO_DIRECTED_ROUNDING", self.payload["status"])

    def test_halo_replay_counts_and_zero_exclusion(self) -> None:
        self.assertEqual(self.theorem["halo_radius"], 2)
        self.assertEqual(self.theorem["total_unique_halo_replay_count"], 14)
        self.assertEqual(self.theorem["total_missing_neighbor_count_due_to_segment_boundaries"], 10)
        self.assertTrue(self.theorem["all_halo_replayed_cells_exclude_zero"])
        self.assertTrue(self.theorem["all_halo_replayed_cells_have_positive_separation"])
        self.assertGreater(Decimal(self.theorem["minimum_halo_decimal_separation"]), Decimal("0"))
        expected_counts = {"zero_projection_amplitude": (8, 4), "stationary_factor": (6, 6)}
        for family in self.theorem["family_halo_replays"]:
            self.assertEqual(
                (family["unique_halo_replay_count"], family["missing_neighbor_count_due_to_segment_boundaries"]),
                expected_counts[family["family"]],
            )
            self.assertTrue(family["all_halo_replayed_cells_exclude_zero"])
            self.assertTrue(family["all_halo_replayed_cells_have_positive_separation"])
            self.assertGreater(Decimal(family["minimum_halo_decimal_separation"]), Decimal("0"))
            for row in family["halo_replay_rows"]:
                self.assertTrue(row["decimal_interval_excludes_zero"])
                self.assertTrue(row["fresh_decimal_separation_positive"])

    def test_chain_counts_and_inherited_p2474(self) -> None:
        self.assertEqual(self.theorem["p2474_total_extremal_witness_rerun_count_inherited"], 28)
        self.assertTrue(self.theorem["p2474_all_fresh_witness_groups_match_stored_inherited"])
        self.assertTrue(self.theorem["p2474_all_fresh_witness_groups_exclude_zero_inherited"])
        self.assertEqual(self.theorem["p2466_tail_count_inherited_from_p2471"], 6361)
        self.assertEqual(self.theorem["p2469_full_opposite_tail_count_inherited"], 45165)
        self.assertEqual(self.theorem["p2470_remaining_non_tail_count_inherited"], 48320)
        self.assertEqual(self.theorem["p2471_p2459_universe_count_inherited"], 99846)
        self.assertEqual(self.theorem["finite_chain_sum_p2466_p2469_p2470"], 99846)
        self.assertTrue(self.theorem["finite_chain_sum_matches_p2471_universe"])
        self.assertEqual(self.theorem["p2471_missing_cells_inherited"], 0)
        self.assertEqual(self.theorem["p2471_extra_cells_inherited"], 0)
        self.assertTrue(self.theorem["p2471_disjoint_complete_inherited"])

    def test_gatekeepers_and_negative_controls(self) -> None:
        for name, value in self.payload["gatekeeper_checks"].items():
            self.assertTrue(value, name)
        self.assertTrue(self.theorem["finite_replay_chain_critical_minimum_halo_replay_audit_exported"])
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
        self.assertIn("one-cell indexing accident", self.theorem["lay_summary"])
        self.assertIn("P2475/S1425", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2475/S1425", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
