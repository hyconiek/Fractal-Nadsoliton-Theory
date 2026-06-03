from __future__ import annotations

import json
import sys
import unittest
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2489_s1439_strict_pointwise_interval_decimal_p2459_boundary_handoff_collar_cumulative_derivative_barrier_certificate import (
    MD,
    OUT,
    build_payload,
    write_markdown,
)


class P2489StrictPointwiseIntervalDecimalP2459BoundaryHandoffCollarCumulativeDerivativeBarrierCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["strict_pointwise_interval_decimal_p2459_boundary_handoff_collar_cumulative_derivative_barrier_certificate"]["theorem_export"]
        cls.barrier = cls.theorem["boundary_handoff_collar_cumulative_derivative_barrier"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2489")
        self.assertEqual(self.payload["stage_id"], "S1439")
        self.assertIn("CUMULATIVE_DERIVATIVE_BARRIER", self.payload["status"])
        self.assertIn("PROOF_COMPRESSION", self.payload["status"])
        self.assertIn("NO_NEW_REPLAY", self.payload["status"])
        self.assertIn("NO_ROOT_WINDOW", self.payload["status"])
        self.assertIn("NO_GLOBAL_ANALYTIC_MONOTONICITY", self.payload["status"])

    def test_barrier_preconditions_and_numeric_export(self) -> None:
        self.assertTrue(self.theorem["all_barrier_preconditions_met"])
        self.assertTrue(all(self.barrier["barrier_preconditions"].values()))
        self.assertTrue(self.theorem["finite_cumulative_lower_barrier_certificate_exported"])
        self.assertEqual(self.theorem["total_barrier_row_count"], 134)
        self.assertEqual(self.theorem["p2487_total_derivative_rows_inherited"], 134)
        self.assertTrue(self.theorem["p2488_checked_collar_lemma_inherited"])
        self.assertTrue(self.theorem["transport_gain_matches_p2488"])
        self.assertGreater(Decimal(self.theorem["minimum_entry_cumulative_lower_barrier"]), Decimal("0"))
        self.assertGreater(Decimal(self.theorem["minimum_exit_cumulative_lower_barrier"]), Decimal("0"))
        self.assertGreater(Decimal(self.theorem["minimum_positive_derivative_gain_on_row"]), Decimal("0"))
        self.assertGreater(Decimal(self.theorem["right_endpoint_cumulative_lower_barrier"]), Decimal(self.barrier["left_anchor_amplitude_lower_bound"]))

    def test_barrier_rows_are_positive_and_cumulative(self) -> None:
        rows = self.barrier["barrier_rows"]
        self.assertEqual(len(rows), 134)
        for prior, current in zip(rows, rows[1:]):
            self.assertEqual(prior["exit_cumulative_lower_barrier"], current["entry_cumulative_lower_barrier"])
        self.assertTrue(all(row["entry_barrier_positive"] for row in rows))
        self.assertTrue(all(row["derivative_gain_positive"] for row in rows))
        self.assertTrue(all(row["exit_barrier_positive"] for row in rows))
        computed_gain = sum(Decimal(row["certified_derivative_gain_on_row"]) for row in rows)
        self.assertEqual(computed_gain, Decimal(self.theorem["total_certified_derivative_gain_over_collar"]))

    def test_proof_compression_and_coverage_budget(self) -> None:
        self.assertTrue(self.theorem["p2489_is_barrier_proof_compression_not_new_replay"])
        self.assertEqual(self.theorem["p2489_new_decimal_replay_row_count"], 0)
        self.assertEqual(self.theorem["p2489_reused_p2487_p2488_row_count_not_a_coverage_count"], 134)
        self.assertEqual(self.theorem["p2489_reused_row_ratio_not_a_p2459_coverage_fraction"], "134/99846")
        self.assertEqual(self.theorem["targeted_p2489_new_p2459_unreplayed_cell_count"], 0)
        self.assertEqual(self.theorem["targeted_p2489_new_p2459_unreplayed_cell_scope_against_p2459_universe"], "0/99846")
        self.assertEqual(self.theorem["p2471_p2459_universe_count_inherited"], 99846)
        self.assertFalse(self.theorem["full_complement_replay_exported_by_this_certificate"])

    def test_gatekeepers_and_negative_controls(self) -> None:
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertFalse(self.theorem["directed_rounding_interval_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["symbolic_root_exclusion_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["root_window_exclusion_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["global_continuum_root_exclusion_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["global_analytic_monotonicity_theorem_exported_by_this_certificate"])
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
        self.assertIn("P2489/S1439", MD.read_text(encoding="utf-8"))
        self.assertIn("P2489/S1439", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2489/S1439", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
