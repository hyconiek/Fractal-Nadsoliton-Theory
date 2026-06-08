from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2584_s1534_apd_full_moment_finite_support_conditional_uniqueness_audit import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2584APDFullMomentFiniteSupportConditionalUniquenessTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["apd_full_moment_finite_support_conditional_uniqueness_audit"]["theorem_export"]
        cls.audit = cls.theorem["apd_full_moment_finite_support_conditional_uniqueness_audit"]

    def test_identity_and_inheritance(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2584")
        self.assertEqual(self.payload["stage_id"], "S1534")
        self.assertIn("FULL_MOMENT_FINITE_SUPPORT_CONDITIONAL", self.payload["status"])
        self.assertEqual(self.theorem["frontier_atom_under_attack"], "strict_dynamical_source_for_A_P_D")
        self.assertTrue(self.theorem["p2416_apd_value_bridge_inherited"])
        self.assertTrue(self.theorem["p2582_low_order_moment_nonuniqueness_inherited"])
        self.assertTrue(self.theorem["p2583_finite_moment_prefix_nonuniqueness_inherited"])

    def test_fixed_support_conditional_uniqueness_but_support_dependence(self) -> None:
        self.assertEqual(self.audit["support_count"], 3)
        self.assertEqual(self.audit["target_count"], 3)
        self.assertEqual(self.audit["full_moment_orders"], [0, 1, 2, 3])
        self.assertTrue(self.audit["all_fixed_supports_have_vandermonde_rank_four"])
        self.assertTrue(self.audit["all_fixed_support_weights_conditionally_unique"])
        self.assertTrue(self.audit["all_recovered_weights_positive"])
        self.assertTrue(self.audit["all_targets_support_dependent"])
        self.assertTrue(self.audit["all_gram_metrics_positive_definite"])
        self.assertTrue(self.audit["all_minimizers_preserve_nodes_and_boundaries"])
        for witness in self.audit["support_witnesses"]:
            self.assertEqual(witness["vandermonde_rank"], 4)
            self.assertLessEqual(witness["max_abs_recovered_weight_error"], 1.0e-10)
            self.assertTrue(witness["full_moments_conditionally_select_weights_on_fixed_support"])
        for target in self.audit["target_rows"]:
            self.assertGreater(target["distinct_middle_probe_values_after_rounding_1e_minus_12"], 1)
            self.assertTrue(target["all_gram_metrics_positive_definite"])
            self.assertTrue(target["all_minimizers_preserve_apd_nodes"])
            self.assertTrue(target["all_minimizers_satisfy_boundary_targets"])
        self.assertTrue(self.theorem["full_moments_select_weights_only_after_support_is_fixed"])
        self.assertTrue(self.theorem["finite_values_boundaries_and_full_support_moments_do_not_select_support"])

    def test_negative_controls_and_docs(self) -> None:
        self.assertFalse(self.theorem["apd_finite_support_source_exported"])
        self.assertFalse(self.theorem["apd_full_moment_law_source_exported"])
        self.assertFalse(self.theorem["apd_support_selection_source_exported"])
        self.assertFalse(self.theorem["apd_positive_measure_source_exported"])
        self.assertFalse(self.theorem["strict_dynamical_source_for_A_P_D_exported"])
        self.assertFalse(self.theorem["bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_certificate"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2584/S1534", MD.read_text(encoding="utf-8"))
        self.assertIn("P2584/S1534", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2584/S1534", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
