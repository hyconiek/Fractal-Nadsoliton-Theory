from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2500_s1450_phase_normalized_third_derivative_symbolic_identity_audit import (
    MD,
    OUT,
    build_payload,
    write_markdown,
)


class P2500PhaseNormalizedThirdDerivativeSymbolicIdentityAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["phase_normalized_third_derivative_symbolic_identity_audit"]["theorem_export"]
        cls.cert = cls.theorem["third_derivative_symbolic_identity_audit"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2500")
        self.assertEqual(self.payload["stage_id"], "S1450")
        self.assertIn("THIRD_DERIVATIVE_SYMBOLIC_IDENTITY", self.payload["status"])
        self.assertIn("NO_FORMAL_DIRECTED_BACKEND", self.payload["status"])
        self.assertIn("NO_GLOBAL_ANALYTIC_UNIQUENESS", self.payload["status"])
        self.assertIn("NO_SOURCE_EXPORT", self.payload["status"])
        self.assertIn("NO_BRIDGE_ATOM", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_symbolic_residuals(self) -> None:
        self.assertEqual(self.theorem["symbolic_backend"], "sympy")
        self.assertTrue(self.theorem["legacy_third_residual_zero"])
        self.assertTrue(self.theorem["strict_third_residual_zero"])
        self.assertTrue(self.theorem["implicit_chain_residual_zero"])
        self.assertTrue(self.theorem["all_symbolic_residuals_zero"])
        self.assertTrue(self.cert["legacy_third_derivative_identity"]["residual_simplifies_to_zero"])
        self.assertTrue(self.cert["strict_third_derivative_identity"]["residual_simplifies_to_zero"])
        self.assertTrue(self.cert["implicit_third_derivative_chain_identity"]["residual_simplifies_to_zero"])

    def test_inherited_p2499_and_gatekeepers(self) -> None:
        self.assertTrue(self.theorem["p2499_local_uniqueness_inherited"])
        self.assertTrue(self.theorem["p2499_third_interval_negative_inherited"])
        self.assertTrue(self.theorem["formula_provenance_supports_p2499_interval_third_derivative"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

    def test_negative_controls(self) -> None:
        self.assertFalse(self.theorem["formal_directed_rounding_backend_exported"])
        self.assertFalse(self.theorem["global_analytic_inflection_uniqueness_theorem_exported"])
        self.assertFalse(self.theorem["curvature_dynamic_source_exported"])
        self.assertFalse(self.theorem["legacy_to_strict_bridge_atom_exported"])
        self.assertFalse(self.theorem["strict_compression_dynamic_source_exported"])
        self.assertFalse(self.theorem["selector_source_theorem_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_certificate"])
        self.assertFalse(self.theorem["role_transfer_licensed_by_this_certificate"])
        self.assertFalse(self.theorem["toe_closure_exported"])

    def test_rg_audit_and_documentation_markers(self) -> None:
        audit = self.payload["rg_non_duplication_audit"]
        self.assertEqual(audit["tool"], "rg")
        self.assertIn("new_packet", audit["patterns"])
        self.assertGreaterEqual(audit["patterns"]["new_packet"]["count"], 0)
        self.assertIn("P2500/S1450", MD.read_text(encoding="utf-8"))
        self.assertIn("P2500/S1450", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2500/S1450", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
