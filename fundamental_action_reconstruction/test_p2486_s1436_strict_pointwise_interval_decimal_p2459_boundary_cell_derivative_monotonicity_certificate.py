from __future__ import annotations

import json
import sys
import unittest
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2486_s1436_strict_pointwise_interval_decimal_p2459_boundary_cell_derivative_monotonicity_certificate import (
    MD,
    OUT,
    build_payload,
    write_markdown,
)


class P2486StrictPointwiseIntervalDecimalP2459BoundaryCellDerivativeMonotonicityCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["strict_pointwise_interval_decimal_p2459_boundary_cell_derivative_monotonicity_certificate"]["theorem_export"]
        cls.certificate = cls.theorem["boundary_cell_derivative_monotonicity_certificate"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2486")
        self.assertEqual(self.payload["stage_id"], "S1436")
        self.assertIn("LOCAL_ONE_CELL", self.payload["status"])
        self.assertIn("NO_ROOT_WINDOW", self.payload["status"])
        self.assertIn("NO_GLOBAL_ANALYTIC_MONOTONICITY", self.payload["status"])
        self.assertIn("NO_COVERAGE_INCREASE", self.payload["status"])

    def test_parent_binding_and_derivative_sign(self) -> None:
        self.assertEqual(self.theorem["parent_collar_side"], "p2456_right_boundary_band")
        self.assertEqual(self.theorem["parent_boundary_band_index"], 0)
        self.assertEqual(self.theorem["parent_dyadic_subcell_index"], 0)
        self.assertEqual(self.theorem["family"], "zero_projection_amplitude")
        self.assertTrue(self.theorem["amplitude_interval_excludes_zero"])
        self.assertGreater(Decimal(self.theorem["amplitude_interval_value"]["lo"]), Decimal("0"))
        self.assertTrue(self.theorem["derivative_interval_excludes_zero"])
        self.assertTrue(self.theorem["derivative_interval_positive_on_entire_cell"])
        self.assertGreater(Decimal(self.theorem["derivative_lower_bound"]), Decimal("0"))
        self.assertGreater(Decimal(self.theorem["derivative_interval_separation_from_zero"]), Decimal("0"))
        self.assertTrue(self.theorem["finite_interval_monotone_increasing_witness"])
        self.assertTrue(self.theorem["left_endpoint_is_interval_minimum_under_derivative_witness"])
        self.assertTrue(self.theorem["zero_exclusion_reinforced_by_positive_value_and_positive_derivative"])

    def test_coverage_budget_and_inherited_counts(self) -> None:
        self.assertEqual(self.theorem["p2485_extended_level_count_inherited"], 64)
        self.assertTrue(self.theorem["p2485_all_secant_margin_drifts_positive_inherited"])
        self.assertEqual(self.theorem["p2471_p2459_universe_count_inherited"], 99846)
        self.assertEqual(self.theorem["p2486_derivative_interval_row_count_not_a_coverage_count"], 1)
        self.assertEqual(self.theorem["p2486_derivative_interval_row_ratio_not_a_p2459_coverage_fraction"], "1/99846")
        self.assertEqual(self.theorem["targeted_p2486_new_p2459_unreplayed_cell_count"], 0)
        self.assertEqual(self.theorem["targeted_p2486_new_p2459_unreplayed_cell_scope_against_p2459_universe"], "0/99846")
        self.assertTrue(self.theorem["p2486_refines_one_inherited_p2456_covered_boundary_chain_cell"])
        self.assertFalse(self.theorem["full_complement_replay_exported_by_this_certificate"])

    def test_gatekeepers_and_negative_controls(self) -> None:
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertTrue(self.theorem["local_one_cell_interval_monotonicity_certificate_exported"])
        self.assertFalse(self.theorem["global_analytic_monotonicity_theorem_exported_by_this_certificate"])
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
        self.assertIn("P2486/S1436", MD.read_text(encoding="utf-8"))
        self.assertIn("P2486/S1436", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2486/S1436", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
