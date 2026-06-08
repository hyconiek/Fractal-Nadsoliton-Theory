from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2590_s1540_apd_finite_even_moment_shell_interval_nonuniqueness_audit import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2590APDFiniteEvenMomentShellIntervalNonuniquenessTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["apd_finite_even_moment_shell_interval_nonuniqueness_audit"]["theorem_export"]
        cls.audit = cls.theorem["apd_finite_even_moment_shell_interval_nonuniqueness_audit"]

    def test_identity_and_inheritance(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2590")
        self.assertEqual(self.payload["stage_id"], "S1540")
        self.assertIn("FINITE_EVEN_MOMENT_SHELL_INTERVAL", self.payload["status"])
        self.assertEqual(self.theorem["frontier_atom_under_attack"], "strict_dynamical_source_for_A_P_D")
        self.assertTrue(self.theorem["p2416_apd_value_bridge_inherited"])
        self.assertTrue(self.theorem["p2584_fixed_support_conditional_uniqueness_inherited"])
        self.assertTrue(self.theorem["p2589_sixth_moment_shell_nonuniqueness_inherited"])

    def test_interval_shell_nonuniqueness(self) -> None:
        self.assertEqual(self.audit["support_count"], 8)
        self.assertEqual(self.audit["target_count"], 3)
        self.assertEqual(self.audit["expected_cardinality"], 10)
        self.assertTrue(self.audit["vieta_power_sum_certificate"]["all_witnesses_share_power_sums"])
        self.assertEqual(self.audit["vieta_power_sum_certificate"]["distinct_product_parameters"], 8)
        self.assertTrue(self.audit["all_supports_have_expected_cardinality"])
        self.assertTrue(self.audit["all_supports_share_endpoints"])
        self.assertTrue(self.audit["all_supports_pair_symmetric_about_center"])
        self.assertTrue(self.audit["all_supports_share_second_fourth_sixth_shell"])
        self.assertTrue(self.audit["all_fixed_support_weights_conditionally_unique"])
        self.assertTrue(self.audit["all_targets_interval_support_nonunique"])
        self.assertTrue(self.audit["all_gram_metrics_positive_definite"])
        self.assertTrue(self.audit["all_minimizers_preserve_nodes_and_boundaries"])
        for witness in self.audit["support_witnesses"]:
            self.assertEqual(witness["support_cardinality"], 10)
            self.assertTrue(witness["is_pair_symmetric_about_center"])
            self.assertTrue(witness["has_expected_internal_second_fourth_sixth_shell"])
            self.assertEqual(witness["vandermonde_rank"], 10)
            self.assertTrue(witness["product_parameter_is_free_interval_coordinate"])
            self.assertTrue(witness["full_moments_conditionally_select_weights_on_fixed_support"])
        for target in self.audit["target_rows"]:
            self.assertGreater(target["distinct_middle_probe_values_after_rounding_1e_minus_12"], 1)
            self.assertTrue(target["all_gram_metrics_positive_definite"])
            self.assertTrue(target["all_minimizers_preserve_apd_nodes"])
            self.assertTrue(target["all_minimizers_satisfy_boundary_targets"])
        self.assertTrue(self.theorem["finite_even_moment_shell_prefix_does_not_select_apd_support"])

    def test_negative_controls_and_docs(self) -> None:
        self.assertFalse(self.theorem["apd_finite_even_moment_shell_source_exported"])
        self.assertFalse(self.theorem["apd_product_parameter_interval_source_exported"])
        self.assertFalse(self.theorem["apd_support_selection_source_exported"])
        self.assertFalse(self.theorem["apd_finite_support_source_exported"])
        self.assertFalse(self.theorem["strict_dynamical_source_for_A_P_D_exported"])
        self.assertFalse(self.theorem["bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_certificate"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2590/S1540", MD.read_text(encoding="utf-8"))
        self.assertIn("P2590/S1540", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2590/S1540", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
