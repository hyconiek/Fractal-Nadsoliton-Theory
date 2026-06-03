from __future__ import annotations

import json
import sys
import unittest
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2488_s1438_strict_pointwise_interval_decimal_p2459_boundary_handoff_collar_monotonicity_lemma_certificate import (
    MD,
    OUT,
    build_payload,
    write_markdown,
)


class P2488StrictPointwiseIntervalDecimalP2459BoundaryHandoffCollarMonotonicityLemmaCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["strict_pointwise_interval_decimal_p2459_boundary_handoff_collar_monotonicity_lemma_certificate"]["theorem_export"]
        cls.lemma = cls.theorem["boundary_handoff_collar_monotonicity_lemma"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2488")
        self.assertEqual(self.payload["stage_id"], "S1438")
        self.assertIn("PROOF_COMPRESSION", self.payload["status"])
        self.assertIn("NO_NEW_REPLAY", self.payload["status"])
        self.assertIn("NO_ROOT_WINDOW", self.payload["status"])
        self.assertIn("NO_GLOBAL_ANALYTIC_MONOTONICITY", self.payload["status"])

    def test_lemma_preconditions_and_exports(self) -> None:
        self.assertTrue(self.theorem["all_lemma_preconditions_met"])
        self.assertTrue(all(self.lemma["lemma_preconditions"].values()))
        self.assertTrue(self.theorem["finite_piecewise_monotone_increasing_collar_lemma_exported"])
        self.assertTrue(self.theorem["finite_positive_collar_zero_exclusion_lemma_exported"])
        self.assertEqual(self.theorem["boundary_band_row_count"], 6)
        self.assertEqual(self.theorem["dyadic_subcell_row_count"], 128)
        self.assertEqual(self.theorem["total_collar_row_count"], 134)
        self.assertGreater(Decimal(self.theorem["minimum_row_amplitude_lower_bound"]), Decimal("0"))
        self.assertGreater(Decimal(self.theorem["minimum_derivative_lower_bound"]), Decimal("0"))
        self.assertGreater(Decimal(self.theorem["minimum_consecutive_positive_amplitude_separation_delta"]), Decimal("0"))

    def test_proof_compression_and_coverage_budget(self) -> None:
        self.assertTrue(self.theorem["p2488_is_proof_compression_not_new_replay"])
        self.assertEqual(self.theorem["p2488_new_decimal_replay_row_count"], 0)
        self.assertEqual(self.theorem["p2488_reused_p2487_derivative_row_count_not_a_coverage_count"], 134)
        self.assertEqual(self.theorem["p2488_reused_derivative_row_ratio_not_a_p2459_coverage_fraction"], "134/99846")
        self.assertEqual(self.theorem["targeted_p2488_new_p2459_unreplayed_cell_count"], 0)
        self.assertEqual(self.theorem["targeted_p2488_new_p2459_unreplayed_cell_scope_against_p2459_universe"], "0/99846")
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
        self.assertIn("P2488/S1438", MD.read_text(encoding="utf-8"))
        self.assertIn("P2488/S1438", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2488/S1438", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
