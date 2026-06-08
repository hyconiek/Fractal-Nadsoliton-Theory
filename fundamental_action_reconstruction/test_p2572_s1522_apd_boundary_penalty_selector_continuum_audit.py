from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2572_s1522_apd_boundary_penalty_selector_continuum_audit import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2572APDBoundaryPenaltySelectorContinuumTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["apd_boundary_penalty_selector_continuum_audit"]["theorem_export"]
        cls.audit = cls.theorem["boundary_penalty_selector_continuum_audit"]

    def test_identity_and_inheritance(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2572")
        self.assertEqual(self.payload["stage_id"], "S1522")
        self.assertIn("BOUNDARY_PENALTY_SELECTOR_CONTINUUM", self.payload["status"])
        self.assertEqual(self.theorem["frontier_atom_under_attack"], "strict_dynamical_source_for_A_P_D")
        self.assertTrue(self.theorem["p2416_apd_value_bridge_inherited"])
        self.assertTrue(self.theorem["p2570_order_obligation_inherited"])
        self.assertTrue(self.theorem["p2571_measure_boundary_obligation_inherited"])

    def test_boundary_penalty_continuum(self) -> None:
        self.assertEqual(len(self.theorem["apd_node_rows"]), 12)
        self.assertEqual(self.audit["fixed_derivative_order"], 2)
        self.assertIn("lambda(t)", self.audit["stationarity_formula"])
        self.assertGreater(self.audit["left_distinct_minimizers_after_rounding_1e_minus_24"], 1)
        self.assertGreater(self.audit["right_distinct_minimizers_after_rounding_1e_minus_24"], 1)
        self.assertTrue(self.audit["left_penalty_changes_selector"])
        self.assertTrue(self.audit["right_penalty_changes_selector"])
        self.assertTrue(self.audit["all_penalty_rows_preserve_apd_nodes"])
        self.assertTrue(self.theorem["finite_apd_values_do_not_select_boundary_penalty"])
        self.assertTrue(self.theorem["fixed_k2_uniform_bulk_still_conditional"])

    def test_negative_controls_and_docs(self) -> None:
        self.assertFalse(self.theorem["apd_boundary_penalty_source_exported"])
        self.assertFalse(self.theorem["apd_boundary_selector_source_exported"])
        self.assertFalse(self.theorem["strict_dynamical_source_for_A_P_D_exported"])
        self.assertFalse(self.theorem["bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_certificate"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2572/S1522", MD.read_text(encoding="utf-8"))
        self.assertIn("P2572/S1522", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2572/S1522", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
