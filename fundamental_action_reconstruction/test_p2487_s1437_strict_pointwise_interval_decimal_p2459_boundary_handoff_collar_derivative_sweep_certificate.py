from __future__ import annotations

import json
import sys
import unittest
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2487_s1437_strict_pointwise_interval_decimal_p2459_boundary_handoff_collar_derivative_sweep_certificate import (
    MD,
    OUT,
    build_payload,
    write_markdown,
)


class P2487StrictPointwiseIntervalDecimalP2459BoundaryHandoffCollarDerivativeSweepCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["strict_pointwise_interval_decimal_p2459_boundary_handoff_collar_derivative_sweep_certificate"]["theorem_export"]
        cls.sweep = cls.theorem["boundary_handoff_collar_derivative_sweep"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2487")
        self.assertEqual(self.payload["stage_id"], "S1437")
        self.assertIn("FINITE_PIECEWISE_MONOTONE", self.payload["status"])
        self.assertIn("NO_ROOT_WINDOW", self.payload["status"])
        self.assertIn("NO_GLOBAL_ANALYTIC_MONOTONICITY", self.payload["status"])
        self.assertIn("NO_COVERAGE_INCREASE", self.payload["status"])

    def test_collar_counts_and_positive_derivatives(self) -> None:
        self.assertEqual(self.theorem["family"], "zero_projection_amplitude")
        self.assertEqual(self.theorem["boundary_band_derivative_row_count"], 6)
        self.assertEqual(self.theorem["dyadic_subcell_derivative_row_count"], 128)
        self.assertEqual(self.theorem["total_derivative_row_count"], 134)
        self.assertEqual(self.theorem["p2481_total_handoff_collar_replay_rows_inherited"], 134)
        self.assertTrue(self.theorem["all_amplitude_intervals_exclude_zero"])
        self.assertTrue(self.theorem["all_amplitude_intervals_positive"])
        self.assertTrue(self.theorem["all_derivative_intervals_exclude_zero"])
        self.assertTrue(self.theorem["all_derivative_intervals_positive"])
        self.assertTrue(self.theorem["all_rows_have_local_interval_monotone_increasing_witness"])
        self.assertGreater(Decimal(self.theorem["minimum_derivative_lower_bound"]), Decimal("0"))

    def test_adjacency_separation_and_piecewise_scope(self) -> None:
        self.assertTrue(self.theorem["all_consecutive_rows_exactly_adjacent"])
        self.assertEqual(Decimal(self.theorem["minimum_endpoint_gap"]), Decimal("0"))
        self.assertEqual(Decimal(self.theorem["maximum_endpoint_gap"]), Decimal("0"))
        self.assertTrue(self.theorem["all_consecutive_amplitude_separations_strictly_increase"])
        self.assertGreater(Decimal(self.theorem["minimum_consecutive_positive_amplitude_separation_delta"]), Decimal("0"))
        self.assertTrue(self.theorem["finite_piecewise_interval_monotonicity_on_checked_collar_exported"])
        self.assertFalse(self.theorem["global_analytic_monotonicity_theorem_exported_by_this_certificate"])

    def test_coverage_budget_and_inherited_counts(self) -> None:
        self.assertTrue(self.theorem["p2486_one_cell_derivative_positive_inherited"])
        self.assertEqual(self.theorem["p2471_p2459_universe_count_inherited"], 99846)
        self.assertEqual(self.theorem["p2487_derivative_interval_row_count_not_a_coverage_count"], 134)
        self.assertEqual(self.theorem["p2487_derivative_interval_row_ratio_not_a_p2459_coverage_fraction"], "134/99846")
        self.assertEqual(self.theorem["targeted_p2487_new_p2459_unreplayed_cell_count"], 0)
        self.assertEqual(self.theorem["targeted_p2487_new_p2459_unreplayed_cell_scope_against_p2459_universe"], "0/99846")
        self.assertTrue(self.theorem["p2487_reuses_p2481_collar_rows_without_new_p2459_coverage"])
        self.assertFalse(self.theorem["full_complement_replay_exported_by_this_certificate"])

    def test_gatekeepers_and_negative_controls(self) -> None:
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertFalse(self.theorem["directed_rounding_interval_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["symbolic_root_exclusion_theorem_exported_by_this_certificate"])
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
        self.assertIn("P2487/S1437", MD.read_text(encoding="utf-8"))
        self.assertIn("P2487/S1437", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2487/S1437", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
