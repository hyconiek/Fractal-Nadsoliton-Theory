from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2503_s1453_strict_damping_marginal_rg_flow_candidate_certificate import (
    MD,
    OUT,
    build_payload,
    write_markdown,
)


class P2503StrictDampingMarginalRGFlowCandidateCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["strict_damping_marginal_rg_flow_candidate_certificate"]["theorem_export"]
        cls.cert = cls.theorem["strict_damping_marginal_rg_flow_candidate_certificate"]
        cls.symbolic = cls.cert["symbolic_rg_flow_candidate"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2503")
        self.assertEqual(self.payload["stage_id"], "S1453")
        self.assertIn("STRICT_DAMPING_MARGINAL_RG_FLOW_CANDIDATE", self.payload["status"])
        self.assertIn("NO_SOURCE_EXPORT", self.payload["status"])
        self.assertIn("NO_BRIDGE_THEOREM", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_symbolic_rg_flow_residuals(self) -> None:
        self.assertEqual(self.theorem["frontier_atom_under_attack"], "strict_damping_beta_eta_source")
        self.assertTrue(self.theorem["frontier_atom_is_in_p2502_bridge_triple"])
        self.assertEqual(self.theorem["gamma_f_plus_delta_minus_eta_residual"], "0")
        self.assertEqual(self.theorem["rg_ode_residual"], "0")
        self.assertEqual(self.theorem["denominator_residual"], "0")
        self.assertTrue(self.theorem["all_symbolic_residuals_zero"])
        self.assertTrue(self.theorem["candidate_solves_formal_denominator_target"])

    def test_delta_omission_and_factorization(self) -> None:
        self.assertTrue(self.theorem["all_delta_omitted_residuals_nonnegative"])
        self.assertGreater(float(self.theorem["max_delta_omitted_residual_on_domain"]), 0.0)
        self.assertLess(float(self.theorem["max_flow_reconstruction_abs_residual_on_domain"]), 1e-70)
        factorization = self.cert["beta_eta_source_candidate_factorization"]
        self.assertIn("delta = 14/5 - 4*log(2)", factorization["required_marginal_correction"])
        self.assertIn("derive delta", factorization["remaining_source_obligation"])

    def test_negative_controls_and_gatekeepers(self) -> None:
        self.assertTrue(self.theorem["candidate_is_not_derived_source_theorem"])
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
        self.assertIn("P2503/S1453", MD.read_text(encoding="utf-8"))
        self.assertIn("P2503/S1453", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2503/S1453", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
