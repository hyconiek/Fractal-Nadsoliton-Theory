from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2582_s1532_apd_low_order_moment_matched_measure_nonuniqueness_audit import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2582APDLowOrderMomentMatchedMeasureNonuniquenessTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["apd_low_order_moment_matched_measure_nonuniqueness_audit"]["theorem_export"]
        cls.audit = cls.theorem["apd_low_order_moment_matched_measure_nonuniqueness_audit"]

    def test_identity_and_inheritance(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2582")
        self.assertEqual(self.payload["stage_id"], "S1532")
        self.assertIn("MOMENT_MATCHED_MEASURE_NONUNIQUENESS", self.payload["status"])
        self.assertEqual(self.theorem["frontier_atom_under_attack"], "strict_dynamical_source_for_A_P_D")
        self.assertTrue(self.theorem["p2416_apd_value_bridge_inherited"])
        self.assertTrue(self.theorem["p2580_metric_covariance_requirement_inherited"])
        self.assertTrue(self.theorem["p2581_measure_dependence_inherited"])

    def test_moment_matched_measure_nonuniqueness(self) -> None:
        self.assertEqual(self.audit["target_count"], 3)
        self.assertEqual(self.audit["moment_matched_measure_count"], 3)
        self.assertEqual(self.audit["moment_orders_matched"], [0, 1, 2])
        self.assertTrue(self.audit["all_measures_positive"])
        self.assertTrue(self.audit["all_measures_share_low_order_moments"])
        self.assertTrue(self.audit["all_targets_low_order_moment_nonunique"])
        self.assertTrue(self.audit["all_gram_metrics_positive_definite"])
        self.assertTrue(self.audit["all_minimizers_preserve_nodes_and_boundaries"])
        for target in self.audit["target_rows"]:
            self.assertEqual(target["moment_matched_measure_count"], 3)
            self.assertGreater(target["distinct_middle_probe_values_after_rounding_1e_minus_12"], 1)
            self.assertTrue(target["all_measures_share_low_order_moments"])
            self.assertTrue(target["all_measures_positive"])
            self.assertTrue(target["all_gram_metrics_positive_definite"])
            self.assertTrue(target["all_minimizers_preserve_apd_nodes"])
            self.assertTrue(target["all_minimizers_satisfy_boundary_targets"])
        self.assertTrue(self.theorem["finite_low_order_moments_do_not_select_positive_measure"])

    def test_negative_controls_and_docs(self) -> None:
        self.assertFalse(self.theorem["apd_low_order_moment_law_source_exported"])
        self.assertFalse(self.theorem["apd_moment_matched_measure_source_exported"])
        self.assertFalse(self.theorem["apd_positive_measure_source_exported"])
        self.assertFalse(self.theorem["apd_l2_inner_product_source_exported"])
        self.assertFalse(self.theorem["strict_dynamical_source_for_A_P_D_exported"])
        self.assertFalse(self.theorem["bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_certificate"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2582/S1532", MD.read_text(encoding="utf-8"))
        self.assertIn("P2582/S1532", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2582/S1532", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
