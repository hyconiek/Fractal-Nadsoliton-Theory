from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2569_s1519_apd_value_bridge_interpolation_dynamic_nonuniqueness_certificate import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2569APDValueBridgeInterpolationDynamicNonuniquenessTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["apd_value_bridge_interpolation_dynamic_nonuniqueness_certificate"]["theorem_export"]
        cls.audit = cls.theorem["interpolation_nonuniqueness_audit"]

    def test_identity_and_inheritance(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2569")
        self.assertEqual(self.payload["stage_id"], "S1519")
        self.assertIn("APD_VALUE_BRIDGE_INTERPOLATION", self.payload["status"])
        self.assertEqual(self.theorem["frontier_atom_under_attack"], "strict_dynamical_source_for_A_P_D")
        self.assertTrue(self.theorem["p2416_apd_value_bridge_inherited"])
        self.assertTrue(self.theorem["p2561_apd_residual_atom_inherited"])
        self.assertTrue(self.theorem["p2568_phase_frequency_still_unsourced_inherited"])

    def test_interpolation_family_nonuniqueness(self) -> None:
        self.assertEqual(len(self.theorem["apd_node_rows"]), 12)
        self.assertEqual(self.audit["base_interpolation_degree"], 11)
        self.assertEqual(self.audit["vanishing_polynomial_degree"], 12)
        self.assertTrue(self.audit["all_family_members_preserve_apd_nodes"])
        self.assertTrue(self.audit["nonzero_family_members_change_off_node_dynamics"])
        self.assertTrue(self.theorem["finite_apd_values_do_not_select_dynamic_law"])
        self.assertTrue(self.theorem["apd_value_bridge_remains_below_dynamic_source"])
        self.assertIn("infinite one-parameter family", self.audit["interpolation_family_cardinality_certified"])

    def test_negative_controls_and_docs(self) -> None:
        self.assertFalse(self.theorem["apd_interpolation_dynamic_source_exported"])
        self.assertFalse(self.theorem["strict_dynamical_source_for_A_P_D_exported"])
        self.assertFalse(self.theorem["bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_certificate"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2569/S1519", MD.read_text(encoding="utf-8"))
        self.assertIn("P2569/S1519", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2569/S1519", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
