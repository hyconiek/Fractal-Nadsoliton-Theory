from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2514_s1464_strict_damping_rg_higher_order_selector_nonidentifiability_certificate import (
    MD,
    OUT,
    build_payload,
    write_markdown,
)


class P2514StrictDampingRGHigherOrderSelectorNonidentifiabilityCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["strict_damping_rg_higher_order_selector_nonidentifiability_certificate"]["theorem_export"]
        cls.cert = cls.theorem["strict_damping_rg_higher_order_selector_nonidentifiability_certificate"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2514")
        self.assertEqual(self.payload["stage_id"], "S1464")
        self.assertIn("HIGHER_ORDER_SELECTOR_NONIDENTIFIABILITY", self.payload["status"])
        self.assertIn("NO_SOURCE_EXPORT", self.payload["status"])
        self.assertIn("NO_BRIDGE_THEOREM", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_higher_order_tower_nonidentifiability(self) -> None:
        self.assertEqual(self.theorem["frontier_atom_under_attack"], "strict_damping_beta_eta_source")
        self.assertTrue(self.theorem["p2513_derivative_order_nonidentifiability_inherited"])
        self.assertEqual(self.theorem["theorem_order_range"], list(range(1, 11)))
        self.assertTrue(self.theorem["all_orders_eliminate_zero_modes_by_node_count"])
        self.assertEqual(self.theorem["finite_audit_orders"], list(range(1, 7)))
        self.assertTrue(self.theorem["finite_gram_audit_supports_orders_1_to_6"])
        self.assertEqual(self.theorem["admissible_order_count_under_node_stationarity_only"], 10)
        self.assertTrue(self.theorem["higher_order_selector_tower_nonidentifiability_exported"])
        finite = self.cert["finite_higher_order_gram_audit"]
        self.assertLess(float(finite["max_abs_node_residual"]), 1e-80)
        for row in finite["order_rows"]:
            self.assertTrue(row["all_leading_minors_positive"])
            self.assertTrue(row["all_pivots_positive"])
            self.assertTrue(row["all_basis_energies_positive"])

    def test_conditional_status_and_negative_controls(self) -> None:
        self.assertTrue(self.theorem["roughness_order_still_requires_strict_source_principle"])
        self.assertTrue(self.theorem["roughness_action_still_postulated_not_derived"])
        self.assertFalse(self.theorem["strict_damping_beta_eta_source_exported"])
        self.assertFalse(self.theorem["strict_dynamical_source_for_A_P_D_exported"])
        self.assertFalse(self.theorem["strict_phase_frequency_source_exported"])
        self.assertFalse(self.theorem["bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["selector_closure_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_certificate"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

    def test_rg_audit_and_documentation_markers(self) -> None:
        audit = self.payload["rg_non_duplication_audit"]
        self.assertEqual(audit["tool"], "rg")
        self.assertIn("new_packet", audit["patterns"])
        self.assertGreaterEqual(audit["patterns"]["new_packet"]["count"], 0)
        self.assertIn("P2514/S1464", MD.read_text(encoding="utf-8"))
        self.assertIn("P2514/S1464", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2514/S1464", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
