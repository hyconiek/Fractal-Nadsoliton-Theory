from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2580_s1530_apd_inner_product_basis_covariance_requirement_certificate import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2580APDInnerProductBasisCovarianceRequirementTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["apd_inner_product_basis_covariance_requirement_certificate"]["theorem_export"]
        cls.audit = cls.theorem["apd_inner_product_basis_covariance_requirement_audit"]

    def test_identity_and_inheritance(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2580")
        self.assertEqual(self.payload["stage_id"], "S1530")
        self.assertIn("BASIS_COVARIANCE_REQUIREMENT", self.payload["status"])
        self.assertEqual(self.theorem["frontier_atom_under_attack"], "strict_dynamical_source_for_A_P_D")
        self.assertTrue(self.theorem["p2416_apd_value_bridge_inherited"])
        self.assertTrue(self.theorem["p2578_basis_metric_dependence_inherited"])
        self.assertTrue(self.theorem["p2579_inverse_metric_tunability_inherited"])

    def test_covariance_requirement(self) -> None:
        self.assertEqual(self.audit["target_count"], 3)
        self.assertEqual(self.audit["basis_variant_count"], 4)
        self.assertTrue(self.audit["naive_euclidean_metric_fails_basis_covariance"])
        self.assertTrue(self.audit["transported_metric_restores_basis_covariance"])
        self.assertTrue(self.audit["transported_metric_is_positive_definite"])
        for target in self.audit["target_rows"]:
            self.assertEqual(target["variant_count"], 4)
            self.assertGreater(target["naive_distinct_middle_probe_values_after_rounding_1e_minus_18"], 1)
            self.assertEqual(target["covariant_distinct_middle_probe_values_after_rounding_1e_minus_6"], 1)
            self.assertTrue(target["naive_euclidean_metric_is_basis_dependent"])
            self.assertTrue(target["covariant_metric_transport_restores_basis_invariance"])
            self.assertTrue(target["covariant_coefficients_match_reference"])
            self.assertTrue(target["all_covariant_metrics_positive_definite"])
        self.assertTrue(self.theorem["coordinate_covariance_removes_basis_artifact_but_not_metric_choice"])

    def test_negative_controls_and_docs(self) -> None:
        self.assertFalse(self.theorem["apd_coordinate_covariant_inner_product_source_exported"])
        self.assertFalse(self.theorem["apd_gram_metric_source_exported"])
        self.assertFalse(self.theorem["apd_basis_covariance_law_source_exported"])
        self.assertFalse(self.theorem["strict_dynamical_source_for_A_P_D_exported"])
        self.assertFalse(self.theorem["bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_certificate"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2580/S1530", MD.read_text(encoding="utf-8"))
        self.assertIn("P2580/S1530", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2580/S1530", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
