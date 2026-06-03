#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
import p2476_s1426_strict_pointwise_interval_decimal_p2459_critical_halo_order_classification_audit as p2476

OUT = ROOT / "generated" / "p2476_s1426_strict_pointwise_interval_decimal_p2459_critical_halo_order_classification_audit.json"
MD = ROOT / "generated" / "p2476_s1426_strict_pointwise_interval_decimal_p2459_critical_halo_order_classification_audit.md"


class P2476StrictPointwiseIntervalDecimalP2459CriticalHaloOrderClassificationAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not OUT.exists():
            p2476.main()
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.theorem = cls.payload["strict_pointwise_interval_decimal_p2459_critical_halo_order_classification_audit"]["theorem_export"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2476")
        self.assertEqual(self.payload["stage_id"], "S1426")
        self.assertIn("P2459_CRITICAL_HALO_ORDER_CLASSIFICATION_AUDIT", self.payload["status"])
        self.assertIn("NO_DIRECTED_ROUNDING", self.payload["status"])

    def test_order_classification_counts(self) -> None:
        self.assertEqual(self.theorem["total_target_classification_count"], 6)
        self.assertEqual(self.theorem["strict_available_halo_minimum_count"], 5)
        self.assertEqual(self.theorem["boundary_truncated_one_sided_minimum_count"], 5)
        self.assertEqual(self.theorem["targets_with_lower_neighbor_count"], 1)
        self.assertFalse(self.theorem["all_targets_are_strict_available_halo_minima"])
        self.assertTrue(self.theorem["lower_neighbor_exception_exported"])
        exceptions = [row for row in self.theorem["target_classifications"] if row["center_has_lower_neighbor_within_available_halo"]]
        self.assertEqual(len(exceptions), 1)
        self.assertEqual(exceptions[0]["family"], "zero_projection_amplitude")
        self.assertEqual(exceptions[0]["source_packet"], "P2469/S1419")
        self.assertGreater(exceptions[0]["lower_neighbor_count"], 0)

    def test_inherited_chain_status(self) -> None:
        self.assertEqual(self.theorem["inherited_p2475_total_unique_halo_replay_count"], 14)
        self.assertTrue(self.theorem["inherited_p2475_all_halos_exclude_zero"])
        self.assertTrue(self.theorem["inherited_p2475_all_halos_positive_separation"])
        self.assertEqual(self.theorem["inherited_p2474_total_extremal_witness_rerun_count"], 28)
        self.assertTrue(self.theorem["inherited_p2473_fingerprint_binding_pass"])
        self.assertTrue(self.theorem["inherited_p2471_disjoint_complete"])

    def test_gatekeepers_and_negative_controls(self) -> None:
        for name, value in self.payload["gatekeeper_checks"].items():
            self.assertTrue(value, name)
        self.assertTrue(self.theorem["finite_replay_chain_critical_halo_order_classification_audit_exported"])
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
        self.assertIn("class-local rather than whole-halo local", self.theorem["lay_summary"])
        self.assertIn("P2476/S1426", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2476/S1426", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
