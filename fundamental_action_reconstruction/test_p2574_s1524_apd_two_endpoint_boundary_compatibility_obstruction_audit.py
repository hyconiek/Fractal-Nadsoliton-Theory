from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2574_s1524_apd_two_endpoint_boundary_compatibility_obstruction_audit import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2574APDTwoEndpointBoundaryCompatibilityObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["apd_two_endpoint_boundary_compatibility_obstruction_audit"]["theorem_export"]
        cls.audit = cls.theorem["two_endpoint_boundary_compatibility_audit"]

    def test_identity_and_inheritance(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2574")
        self.assertEqual(self.payload["stage_id"], "S1524")
        self.assertIn("TWO_ENDPOINT_BOUNDARY_COMPATIBILITY", self.payload["status"])
        self.assertEqual(self.theorem["frontier_atom_under_attack"], "strict_dynamical_source_for_A_P_D")
        self.assertTrue(self.theorem["p2416_apd_value_bridge_inherited"])
        self.assertTrue(self.theorem["p2572_boundary_penalty_continuum_inherited"])
        self.assertTrue(self.theorem["p2573_inverse_boundary_tunability_inherited"])

    def test_two_endpoint_boundary_obstruction(self) -> None:
        self.assertEqual(len(self.theorem["apd_node_rows"]), 12)
        self.assertEqual(self.audit["target_count"], 3)
        self.assertIn("q_lambda'", self.audit["boundary_slope_formula"])
        self.assertTrue(self.audit["all_targets_fail_exact_two_endpoint_compatibility"])
        self.assertTrue(self.audit["zero_zero_neumann_incompatible"])
        self.assertTrue(self.audit["all_least_squares_members_preserve_apd_nodes"])
        zero_row = next(row for row in self.audit["target_rows"] if row["target_name"] == "zero_zero_neumann")
        self.assertGreater(zero_row["endpoint_lambda_gap"], 0.0)
        self.assertTrue(self.theorem["finite_apd_values_do_not_supply_compatible_two_endpoint_boundary_law"])

    def test_negative_controls_and_docs(self) -> None:
        self.assertFalse(self.theorem["apd_two_endpoint_boundary_source_exported"])
        self.assertFalse(self.theorem["apd_neumann_boundary_source_exported"])
        self.assertFalse(self.theorem["strict_dynamical_source_for_A_P_D_exported"])
        self.assertFalse(self.theorem["bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_certificate"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2574/S1524", MD.read_text(encoding="utf-8"))
        self.assertIn("P2574/S1524", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2574/S1524", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
