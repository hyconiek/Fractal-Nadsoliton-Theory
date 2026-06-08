from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2579_s1529_apd_inner_product_inverse_metric_tunability_audit import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2579APDInnerProductInverseMetricTunabilityTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["apd_inner_product_inverse_metric_tunability_audit"]["theorem_export"]
        cls.audit = cls.theorem["apd_inner_product_inverse_metric_tunability_audit"]

    def test_identity_and_inheritance(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2579")
        self.assertEqual(self.payload["stage_id"], "S1529")
        self.assertIn("INVERSE_METRIC_TUNABILITY", self.payload["status"])
        self.assertEqual(self.theorem["frontier_atom_under_attack"], "strict_dynamical_source_for_A_P_D")
        self.assertTrue(self.theorem["p2416_apd_value_bridge_inherited"])
        self.assertTrue(self.theorem["p2575_augmented_nullspace_inherited"])
        self.assertTrue(self.theorem["p2578_basis_metric_dependence_inherited"])

    def test_inverse_metric_tunability(self) -> None:
        self.assertEqual(self.audit["target_count"], 3)
        self.assertEqual(len(self.audit["target_gamma_values"]), 5)
        self.assertEqual(self.audit["boundary_matrix_rank"], 2)
        self.assertEqual(self.audit["coefficient_nullity"], 1)
        self.assertTrue(self.audit["all_targets_inverse_metric_tunable"])
        self.assertTrue(self.audit["all_constructed_metrics_positive_definite"])
        self.assertTrue(self.audit["all_constructed_metrics_preserve_nodes_boundaries_and_stationarity"])
        for target in self.audit["target_rows"]:
            self.assertEqual(target["metric_witness_count"], 5)
            self.assertGreater(target["distinct_middle_probe_values_after_rounding_1e_minus_18"], 1)
            self.assertTrue(target["all_witness_metrics_positive_definite"])
            self.assertTrue(target["all_witnesses_hit_requested_gamma"])
            self.assertTrue(target["all_witnesses_stationary_for_their_metric"])
            self.assertTrue(target["all_witnesses_preserve_apd_nodes"])
            self.assertTrue(target["all_witnesses_satisfy_boundary_targets"])
        self.assertTrue(self.theorem["finite_apd_values_and_boundary_targets_do_not_select_function_space_inner_product"])

    def test_negative_controls_and_docs(self) -> None:
        self.assertFalse(self.theorem["apd_function_space_inner_product_source_exported"])
        self.assertFalse(self.theorem["apd_spd_metric_source_exported"])
        self.assertFalse(self.theorem["apd_inverse_metric_selector_source_exported"])
        self.assertFalse(self.theorem["strict_dynamical_source_for_A_P_D_exported"])
        self.assertFalse(self.theorem["bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_certificate"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2579/S1529", MD.read_text(encoding="utf-8"))
        self.assertIn("P2579/S1529", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2579/S1529", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
