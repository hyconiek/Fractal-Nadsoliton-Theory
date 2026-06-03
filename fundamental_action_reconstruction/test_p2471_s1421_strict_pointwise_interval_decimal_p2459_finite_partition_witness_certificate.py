#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
import unittest
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
import p2471_s1421_strict_pointwise_interval_decimal_p2459_finite_partition_witness_certificate as p2471

OUT = ROOT / "generated" / "p2471_s1421_strict_pointwise_interval_decimal_p2459_finite_partition_witness_certificate.json"
MD = ROOT / "generated" / "p2471_s1421_strict_pointwise_interval_decimal_p2459_finite_partition_witness_certificate.md"


class P2471StrictPointwiseIntervalDecimalP2459FinitePartitionWitnessCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not OUT.exists():
            p2471.main()
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.theorem = cls.payload["strict_pointwise_interval_decimal_p2459_finite_partition_witness_certificate"]["theorem_export"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2471")
        self.assertEqual(self.payload["stage_id"], "S1421")
        self.assertIn("P2459_FINITE_PARTITION_WITNESS", self.payload["status"])
        self.assertIn("NO_DIRECTED_ROUNDING", self.payload["status"])

    def test_inherited_and_partition_counts(self) -> None:
        self.assertEqual(self.theorem["p2459_total_unreplayed_by_decimal_boundary_chain_cell_count_inherited"], 99846)
        self.assertEqual(self.theorem["p2466_total_tail_boundary_replay_count_inherited"], 6361)
        self.assertEqual(self.theorem["p2469_total_opposite_tail_full_replay_count_inherited"], 45165)
        self.assertEqual(self.theorem["p2470_total_remaining_non_tail_full_replay_count_inherited"], 48320)
        self.assertEqual(self.theorem["total_p2459_universe_cell_count_rebuilt"], 99846)
        self.assertEqual(self.theorem["total_p2466_descent_tail_cell_count_in_partition"], 6361)
        self.assertEqual(self.theorem["total_p2469_opposite_tail_cell_count_in_partition"], 45165)
        self.assertEqual(self.theorem["total_p2470_remaining_non_tail_cell_count_in_partition"], 48320)
        self.assertEqual(self.theorem["total_partition_union_cell_count"], 99846)
        self.assertEqual(self.theorem["total_partition_missing_cell_count"], 0)
        self.assertEqual(self.theorem["total_partition_extra_cell_count"], 0)

    def test_family_partitions(self) -> None:
        expected = {
            "zero_projection_amplitude": (49898, 1136, 8269, 40493),
            "stationary_factor": (49948, 5225, 36896, 7827),
        }
        for row in self.theorem["family_partition_witnesses"]:
            self.assertTrue(row["is_disjoint_partition_of_p2459_unreplayed_by_boundary_chain_universe"])
            self.assertEqual(row["partition_missing_cell_count"], 0)
            self.assertEqual(row["partition_extra_cell_count"], 0)
            self.assertTrue(all(value == 0 for value in row["pairwise_intersection_counts"].values()))
            self.assertTrue(all(row["inherited_counts_match_independent_partition_sets"].values()))
            self.assertGreater(Decimal(row["minimum_inherited_decimal_separation_across_partition"]), Decimal("0"))
            self.assertEqual(
                (
                    row["p2459_universe_cell_count"],
                    row["p2466_descent_tail_cell_count"],
                    row["p2469_opposite_tail_cell_count"],
                    row["p2470_remaining_non_tail_cell_count"],
                ),
                expected[row["family"]],
            )

    def test_gatekeepers_and_negative_controls(self) -> None:
        for name, value in self.payload["gatekeeper_checks"].items():
            self.assertTrue(value, name)
        self.assertTrue(self.theorem["finite_p2459_partition_witness_exported"])
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

    def test_rg_audit_docs_and_lay_summary(self) -> None:
        audit = self.payload["rg_non_duplication_audit"]
        self.assertEqual(audit["tool"], "rg")
        self.assertIn("new_packet", audit["patterns"])
        md_text = MD.read_text(encoding="utf-8")
        self.assertIn("Plain-language progress note", md_text)
        self.assertIn("no box is missing", self.theorem["lay_summary"])
        self.assertIn("P2471/S1421", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2471/S1421", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
