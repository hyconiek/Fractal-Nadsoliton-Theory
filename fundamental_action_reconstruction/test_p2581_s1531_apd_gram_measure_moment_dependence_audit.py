from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2581_s1531_apd_gram_measure_moment_dependence_audit import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2581APDGramMeasureMomentDependenceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["apd_gram_measure_moment_dependence_audit"]["theorem_export"]
        cls.audit = cls.theorem["apd_gram_measure_moment_dependence_audit"]

    def test_identity_and_inheritance(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2581")
        self.assertEqual(self.payload["stage_id"], "S1531")
        self.assertIn("GRAM_MEASURE_MOMENT_DEPENDENCE", self.payload["status"])
        self.assertEqual(self.theorem["frontier_atom_under_attack"], "strict_dynamical_source_for_A_P_D")
        self.assertTrue(self.theorem["p2416_apd_value_bridge_inherited"])
        self.assertTrue(self.theorem["p2579_inner_product_tunability_inherited"])
        self.assertTrue(self.theorem["p2580_metric_covariance_requirement_inherited"])

    def test_measure_moment_dependence(self) -> None:
        self.assertEqual(self.audit["target_count"], 3)
        self.assertEqual(self.audit["measure_variant_count"], 5)
        self.assertTrue(self.audit["all_targets_measure_dependent"])
        self.assertTrue(self.audit["all_gram_metrics_positive_definite"])
        self.assertTrue(self.audit["all_measure_minimizers_preserve_nodes_and_boundaries"])
        for target in self.audit["target_rows"]:
            self.assertEqual(target["measure_count"], 5)
            self.assertGreater(target["distinct_middle_probe_values_after_rounding_1e_minus_12"], 1)
            self.assertTrue(target["all_gram_metrics_positive_definite"])
            self.assertTrue(target["all_measure_minimizers_preserve_apd_nodes"])
            self.assertTrue(target["all_measure_minimizers_satisfy_boundary_targets"])
        self.assertTrue(self.theorem["finite_apd_values_and_boundary_targets_do_not_select_positive_measure"])

    def test_negative_controls_and_docs(self) -> None:
        self.assertFalse(self.theorem["apd_positive_measure_source_exported"])
        self.assertFalse(self.theorem["apd_gram_moment_source_exported"])
        self.assertFalse(self.theorem["apd_l2_inner_product_source_exported"])
        self.assertFalse(self.theorem["apd_measure_selector_source_exported"])
        self.assertFalse(self.theorem["strict_dynamical_source_for_A_P_D_exported"])
        self.assertFalse(self.theorem["bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_certificate"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2581/S1531", MD.read_text(encoding="utf-8"))
        self.assertIn("P2581/S1531", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2581/S1531", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
