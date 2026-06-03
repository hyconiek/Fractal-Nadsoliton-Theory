from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2510_s1460_strict_damping_rg_roughness_kkt_stationarity_certificate import (
    MD,
    OUT,
    build_payload,
    write_markdown,
)


class P2510StrictDampingRGRoughnessKKTStationarityCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["strict_damping_rg_roughness_kkt_stationarity_certificate"]["theorem_export"]
        cls.cert = cls.theorem["strict_damping_rg_roughness_kkt_stationarity_certificate"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2510")
        self.assertEqual(self.payload["stage_id"], "S1460")
        self.assertIn("ROUGHNESS_KKT_STATIONARITY", self.payload["status"])
        self.assertIn("NO_SOURCE_EXPORT", self.payload["status"])
        self.assertIn("NO_BRIDGE_THEOREM", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_symbolic_and_finite_kkt_stationarity(self) -> None:
        self.assertEqual(self.theorem["frontier_atom_under_attack"], "strict_damping_beta_eta_source")
        self.assertTrue(self.theorem["p2509_variational_wellposedness_inherited"])
        self.assertTrue(self.theorem["symbolic_stationarity_residual_zero"])
        self.assertTrue(self.theorem["finite_kkt_all_full_rank"])
        self.assertTrue(self.theorem["finite_kkt_all_affine_zero_multipliers"])
        self.assertTrue(self.theorem["finite_kkt_all_constraints_satisfied"])
        self.assertTrue(self.theorem["finite_kkt_all_stationarity_residuals_zero"])
        self.assertTrue(self.theorem["kkt_stationarity_confirmed_for_postulated_functional"])
        rows = self.cert["finite_polynomial_kkt_audit"]["degree_rows"]
        self.assertEqual([row["polynomial_degree"] for row in rows], [12, 14])
        for row in rows:
            self.assertTrue(row["kkt_full_rank"])
            self.assertTrue(row["solution_is_affine_delta_ell_with_zero_multipliers"])
            self.assertLess(float(row["max_abs_linear_solve_residual"]), 1e-80)
            self.assertLess(float(row["max_abs_coeff_error_vs_affine_delta_ell"]), 1e-70)
            self.assertLess(float(row["max_abs_node_multiplier"]), 1e-70)

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
        self.assertIn("P2510/S1460", MD.read_text(encoding="utf-8"))
        self.assertIn("P2510/S1460", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2510/S1460", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
