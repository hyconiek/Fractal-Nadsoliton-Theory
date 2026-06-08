from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2577_s1527_apd_boundary_nullspace_grid_measure_dependence_audit import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2577APDBoundaryNullspaceGridMeasureDependenceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["apd_boundary_nullspace_grid_measure_dependence_audit"]["theorem_export"]
        cls.audit = cls.theorem["boundary_nullspace_grid_measure_dependence_audit"]

    def test_identity_and_inheritance(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2577")
        self.assertEqual(self.payload["stage_id"], "S1527")
        self.assertIn("BOUNDARY_NULLSPACE_GRID_MEASURE_DEPENDENCE", self.payload["status"])
        self.assertEqual(self.theorem["frontier_atom_under_attack"], "strict_dynamical_source_for_A_P_D")
        self.assertTrue(self.theorem["p2416_apd_value_bridge_inherited"])
        self.assertTrue(self.theorem["p2575_augmented_nullspace_inherited"])
        self.assertTrue(self.theorem["p2576_order_dependence_inherited"])

    def test_grid_measure_dependence(self) -> None:
        self.assertEqual(self.audit["boundary_matrix_rank"], 2)
        self.assertEqual(self.audit["boundary_matrix_nullity"], 1)
        self.assertEqual(self.audit["target_count"], 3)
        self.assertEqual(self.audit["variant_count_per_target"], 7)
        self.assertEqual(self.audit["fixed_derivative_order"], 2)
        self.assertTrue(self.audit["all_targets_have_grid_measure_dependent_gamma_selector"])
        self.assertTrue(self.audit["all_selected_gammas_preserve_nodes_and_boundaries"])
        for target in self.audit["target_rows"]:
            self.assertGreater(target["distinct_gamma_minimizers_after_rounding_1e_minus_20"], 1)
            self.assertTrue(target["all_gamma_minimizers_preserve_apd_nodes"])
            self.assertTrue(target["all_gamma_minimizers_preserve_boundary_targets"])
        self.assertTrue(self.theorem["finite_apd_values_and_boundary_targets_do_not_select_grid_measure"])

    def test_negative_controls_and_docs(self) -> None:
        self.assertFalse(self.theorem["apd_boundary_nullspace_grid_measure_source_exported"])
        self.assertFalse(self.theorem["apd_nullspace_gamma_selector_source_exported"])
        self.assertFalse(self.theorem["apd_quadrature_measure_source_exported"])
        self.assertFalse(self.theorem["strict_dynamical_source_for_A_P_D_exported"])
        self.assertFalse(self.theorem["bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_certificate"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2577/S1527", MD.read_text(encoding="utf-8"))
        self.assertIn("P2577/S1527", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2577/S1527", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
