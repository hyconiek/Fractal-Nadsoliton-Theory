from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2575_s1525_apd_augmented_boundary_nullspace_nonuniqueness_audit import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2575APDAugmentedBoundaryNullspaceNonuniquenessTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["apd_augmented_boundary_nullspace_nonuniqueness_audit"]["theorem_export"]
        cls.audit = cls.theorem["augmented_boundary_nullspace_audit"]

    def test_identity_and_inheritance(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2575")
        self.assertEqual(self.payload["stage_id"], "S1525")
        self.assertIn("AUGMENTED_BOUNDARY_NULLSPACE_NONUNIQUENESS", self.payload["status"])
        self.assertEqual(self.theorem["frontier_atom_under_attack"], "strict_dynamical_source_for_A_P_D")
        self.assertTrue(self.theorem["p2416_apd_value_bridge_inherited"])
        self.assertTrue(self.theorem["p2573_inverse_boundary_tunability_inherited"])
        self.assertTrue(self.theorem["p2574_two_endpoint_obstruction_inherited"])

    def test_augmented_boundary_nullspace(self) -> None:
        self.assertEqual(len(self.theorem["apd_node_rows"]), 12)
        self.assertEqual(self.audit["boundary_matrix_rank"], 2)
        self.assertEqual(self.audit["boundary_matrix_nullity"], 1)
        self.assertEqual(self.audit["target_count"], 3)
        self.assertTrue(self.audit["all_targets_solved_by_augmented_basis"])
        self.assertTrue(self.audit["all_targets_preserve_finite_apd_nodes"])
        self.assertTrue(self.audit["nullspace_keeps_boundary_targets_but_changes_off_node_values"])
        self.assertTrue(self.theorem["two_endpoint_boundary_can_be_satisfied_after_adding_vanishing_modes"])
        self.assertTrue(self.theorem["but_augmented_boundary_solution_is_nonunique"])

    def test_negative_controls_and_docs(self) -> None:
        self.assertFalse(self.theorem["apd_augmented_boundary_source_exported"])
        self.assertFalse(self.theorem["apd_vanishing_mode_source_exported"])
        self.assertFalse(self.theorem["apd_boundary_nullspace_selector_source_exported"])
        self.assertFalse(self.theorem["strict_dynamical_source_for_A_P_D_exported"])
        self.assertFalse(self.theorem["bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_certificate"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2575/S1525", MD.read_text(encoding="utf-8"))
        self.assertIn("P2575/S1525", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2575/S1525", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
