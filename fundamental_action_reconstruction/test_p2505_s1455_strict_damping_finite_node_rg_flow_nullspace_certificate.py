from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2505_s1455_strict_damping_finite_node_rg_flow_nullspace_certificate import (
    MD,
    OUT,
    build_payload,
    write_markdown,
)


class P2505StrictDampingFiniteNodeRGFlowNullspaceCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["strict_damping_finite_node_rg_flow_nullspace_certificate"]["theorem_export"]
        cls.cert = cls.theorem["strict_damping_finite_node_rg_flow_nullspace_certificate"]
        cls.symbolic = cls.cert["symbolic_nullspace_certificate"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2505")
        self.assertEqual(self.payload["stage_id"], "S1455")
        self.assertIn("FINITE_NODE_RG_FLOW_NULLSPACE", self.payload["status"])
        self.assertIn("NO_SOURCE_EXPORT", self.payload["status"])
        self.assertIn("NO_BRIDGE_THEOREM", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_symbolic_nullspace(self) -> None:
        self.assertEqual(self.theorem["frontier_atom_under_attack"], "strict_damping_beta_eta_source")
        self.assertTrue(self.theorem["p2504_constant_coefficient_uniqueness_inherited"])
        self.assertTrue(self.theorem["all_node_residuals_zero_symbolically"])
        self.assertEqual(len(self.symbolic["node_rows"]), 11)
        self.assertTrue(all(row["R_log_d"] == "0" for row in self.symbolic["node_rows"]))
        self.assertTrue(self.theorem["flow_lambda_nonconstant_for_nonzero_epsilon"])

    def test_finite_node_and_midpoint_witness(self) -> None:
        self.assertLess(float(self.theorem["max_abs_node_reconstruction_residual"]), 1e-70)
        self.assertGreater(float(self.theorem["max_abs_midpoint_lambda_deviation"]), 0.0)
        self.assertTrue(self.theorem["all_finite_nodes_preserved_within_precision"])
        self.assertTrue(self.theorem["midpoint_flow_deviation_nonzero"])
        self.assertTrue(self.theorem["finite_nodes_do_not_select_constant_coefficient_flow"])
        self.assertTrue(self.theorem["p2504_uniqueness_limited_to_constant_coefficient_ansatz_confirmed"])

    def test_negative_controls_and_gatekeepers(self) -> None:
        self.assertTrue(self.theorem["nonconstant_nullspace_not_source_theorem"])
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
        self.assertIn("P2505/S1455", MD.read_text(encoding="utf-8"))
        self.assertIn("P2505/S1455", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2505/S1455", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
