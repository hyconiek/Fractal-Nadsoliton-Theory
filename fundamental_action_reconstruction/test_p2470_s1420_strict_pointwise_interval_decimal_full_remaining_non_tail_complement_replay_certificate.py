#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
import unittest
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
import p2470_s1420_strict_pointwise_interval_decimal_full_remaining_non_tail_complement_replay_certificate as p2470

OUT = ROOT / "generated" / "p2470_s1420_strict_pointwise_interval_decimal_full_remaining_non_tail_complement_replay_certificate.json"
MD = ROOT / "generated" / "p2470_s1420_strict_pointwise_interval_decimal_full_remaining_non_tail_complement_replay_certificate.md"


class P2470StrictPointwiseIntervalDecimalFullRemainingNonTailComplementReplayCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not OUT.exists():
            p2470.main()
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.theorem = cls.payload["strict_pointwise_interval_decimal_full_remaining_non_tail_complement_replay_certificate"]["theorem_export"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2470")
        self.assertEqual(self.payload["stage_id"], "S1420")
        self.assertIn("FULL_REMAINING_NON_TAIL_COMPLEMENT_REPLAY", self.payload["status"])
        self.assertIn("FINITE_GRID_ONLY", self.payload["status"])

    def test_inherited_counts(self) -> None:
        self.assertEqual(self.theorem["p2451_total_interval_complement_cell_count_inherited"], 99882)
        self.assertEqual(self.theorem["p2459_total_unreplayed_by_decimal_boundary_chain_cell_count_inherited"], 99846)
        self.assertEqual(self.theorem["p2466_total_tail_boundary_replay_count_inherited"], 6361)
        self.assertEqual(self.theorem["p2469_total_opposite_tail_full_replay_count_inherited"], 45165)
        self.assertEqual(self.theorem["p2469_p2459_unreplayed_budget_not_covered_by_tail_replays_inherited"], 48320)

    def test_remaining_non_tail_and_p2459_accounting(self) -> None:
        self.assertEqual(self.theorem["total_remaining_non_tail_available_cell_count"], 48320)
        self.assertEqual(self.theorem["total_remaining_non_tail_full_replay_count"], 48320)
        self.assertEqual(self.theorem["total_remaining_non_tail_unreplayed_budget_after_full_replay"], 0)
        self.assertEqual(self.theorem["p2459_total_finite_decimal_replayed_cell_count_from_p2466_p2469_p2470"], 99846)
        self.assertEqual(self.theorem["p2459_unreplayed_by_boundary_chain_budget_after_p2466_p2469_p2470"], 0)
        self.assertEqual(Decimal(self.theorem["p2459_unreplayed_by_boundary_chain_finite_decimal_coverage_ratio"]), Decimal("1"))
        self.assertTrue(self.theorem["finite_decimal_replay_covers_full_p2459_unreplayed_by_boundary_chain_budget"])

    def test_zero_exclusion_and_minimum(self) -> None:
        self.assertTrue(self.theorem["all_remaining_non_tail_full_replayed_cells_exclude_zero"])
        self.assertTrue(self.theorem["all_remaining_non_tail_replay_keys_unique_by_family"])
        self.assertGreater(Decimal(self.theorem["minimum_remaining_non_tail_full_replay_decimal_separation"]), Decimal("0"))
        family_counts = {row["family"]: row["remaining_non_tail_full_replay_count"] for row in self.theorem["family_remaining_non_tail_replays"]}
        self.assertEqual(family_counts["zero_projection_amplitude"], 40493)
        self.assertEqual(family_counts["stationary_factor"], 7827)
        for family in self.theorem["family_remaining_non_tail_replays"]:
            self.assertEqual(family["remaining_non_tail_unreplayed_budget_after_full_replay"], 0)
            self.assertEqual(family["p2459_family_budget_after_p2466_p2469_p2470_replays"], 0)
            self.assertTrue(family["all_remaining_non_tail_full_replayed_cells_exclude_zero"])
            self.assertGreater(Decimal(family["minimum_remaining_non_tail_full_replay_decimal_separation"]), Decimal("0"))

    def test_gatekeepers_and_negative_controls(self) -> None:
        for name, value in self.payload["gatekeeper_checks"].items():
            self.assertTrue(value, name)
        self.assertTrue(self.theorem["finite_grid_complement_ledger_exported"])
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

    def test_rg_audit_docs_and_lay_note(self) -> None:
        audit = self.payload["rg_non_duplication_audit"]
        self.assertEqual(audit["tool"], "rg")
        self.assertIn("new_packet", audit["patterns"])
        md_text = MD.read_text(encoding="utf-8")
        self.assertIn("Plain-language progress note", md_text)
        self.assertIn("P2470/S1420", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2470/S1420", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
