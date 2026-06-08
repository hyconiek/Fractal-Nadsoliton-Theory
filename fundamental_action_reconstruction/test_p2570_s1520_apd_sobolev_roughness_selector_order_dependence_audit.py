from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2570_s1520_apd_sobolev_roughness_selector_order_dependence_audit import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2570APDSobolevRoughnessSelectorOrderDependenceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["apd_sobolev_roughness_selector_order_dependence_audit"]["theorem_export"]
        cls.audit = cls.theorem["sobolev_roughness_order_dependence_audit"]

    def test_identity_and_inheritance(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2570")
        self.assertEqual(self.payload["stage_id"], "S1520")
        self.assertIn("SOBOLEV_ROUGHNESS_SELECTOR_ORDER_DEPENDENCE", self.payload["status"])
        self.assertEqual(self.theorem["frontier_atom_under_attack"], "strict_dynamical_source_for_A_P_D")
        self.assertTrue(self.theorem["p2416_apd_value_bridge_inherited"])
        self.assertTrue(self.theorem["p2561_apd_residual_atom_inherited"])
        self.assertTrue(self.theorem["p2569_interpolation_family_inherited"])

    def test_roughness_order_dependence(self) -> None:
        self.assertEqual(len(self.theorem["apd_node_rows"]), 12)
        self.assertEqual(self.audit["derivative_orders"], [0, 1, 2, 3, 4, 12])
        self.assertEqual(len(self.audit["roughness_audit_rows"]), 6)
        self.assertTrue(self.audit["all_selected_members_preserve_apd_nodes"])
        self.assertTrue(self.audit["roughness_order_changes_selector"])
        self.assertIn(12, self.audit["orders_selecting_base_lambda_zero"])
        self.assertIn(0, self.audit["orders_selecting_nonzero_lambda"])
        self.assertGreater(self.audit["distinct_minimizers_after_rounding_1e_minus_24"], 1)
        self.assertTrue(self.theorem["finite_apd_values_do_not_select_sobolev_order"])
        self.assertTrue(self.theorem["sobolev_selector_is_conditional_not_strict_source"])

    def test_negative_controls_and_docs(self) -> None:
        self.assertFalse(self.theorem["apd_sobolev_order_source_exported"])
        self.assertFalse(self.theorem["apd_roughness_selector_source_exported"])
        self.assertFalse(self.theorem["strict_dynamical_source_for_A_P_D_exported"])
        self.assertFalse(self.theorem["bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_certificate"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2570/S1520", MD.read_text(encoding="utf-8"))
        self.assertIn("P2570/S1520", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2570/S1520", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
