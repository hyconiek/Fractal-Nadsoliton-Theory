from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2513_s1463_strict_damping_rg_derivative_order_nonidentifiability_certificate import (
    MD,
    OUT,
    build_payload,
    write_markdown,
)


class P2513StrictDampingRGDerivativeOrderNonidentifiabilityCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["strict_damping_rg_derivative_order_nonidentifiability_certificate"]["theorem_export"]
        cls.cert = cls.theorem["strict_damping_rg_derivative_order_nonidentifiability_certificate"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2513")
        self.assertEqual(self.payload["stage_id"], "S1463")
        self.assertIn("DERIVATIVE_ORDER_NONIDENTIFIABILITY", self.payload["status"])
        self.assertIn("NO_SOURCE_EXPORT", self.payload["status"])
        self.assertIn("NO_BRIDGE_THEOREM", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_derivative_order_nonidentifiability(self) -> None:
        self.assertEqual(self.theorem["frontier_atom_under_attack"], "strict_damping_beta_eta_source")
        self.assertTrue(self.theorem["p2512_derivative_only_ambiguity_inherited"])
        self.assertTrue(self.theorem["h1_and_h2_both_select_affine_node_solution"])
        self.assertTrue(self.theorem["h1_gram_positive_on_finite_tangent_basis"])
        self.assertTrue(self.theorem["h2_gram_positive_on_finite_tangent_basis"])
        self.assertTrue(self.theorem["all_mixed_nonnegative_nonzero_rows_positive"])
        self.assertTrue(self.theorem["finite_gram_supports_h1_h2_and_mixed_coercivity"])
        self.assertTrue(self.theorem["derivative_order_nonidentifiability_exported"])
        self.assertTrue(self.theorem["roughness_order_still_requires_strict_source_principle"])
        finite = self.cert["finite_derivative_order_audit"]
        self.assertLess(float(finite["max_abs_node_residual"]), 1e-80)
        self.assertTrue(finite["all_basis_rows_positive_for_h1_and_h2"])

    def test_conditional_status_and_negative_controls(self) -> None:
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
        self.assertIn("P2513/S1463", MD.read_text(encoding="utf-8"))
        self.assertIn("P2513/S1463", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2513/S1463", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
