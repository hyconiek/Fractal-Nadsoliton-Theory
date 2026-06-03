from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2504_s1454_strict_damping_constant_coefficient_rg_uniqueness_certificate import (
    MD,
    OUT,
    build_payload,
    write_markdown,
)


class P2504StrictDampingConstantCoefficientRGUniquenessCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["strict_damping_constant_coefficient_rg_uniqueness_certificate"]["theorem_export"]
        cls.cert = cls.theorem["strict_damping_constant_coefficient_rg_uniqueness_certificate"]
        cls.symbolic = cls.cert["symbolic_two_node_uniqueness"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2504")
        self.assertEqual(self.payload["stage_id"], "S1454")
        self.assertIn("CONSTANT_COEFFICIENT_RG_UNIQUENESS", self.payload["status"])
        self.assertIn("NO_SOURCE_EXPORT", self.payload["status"])
        self.assertIn("NO_BRIDGE_THEOREM", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_symbolic_two_node_solution(self) -> None:
        self.assertEqual(self.theorem["frontier_atom_under_attack"], "strict_damping_beta_eta_source")
        self.assertTrue(self.theorem["p2503_candidate_inherited"])
        self.assertEqual(self.theorem["beta0_solution"], "1")
        self.assertEqual(self.theorem["lambda_solution"], "14/5 - 4*log(2)")
        self.assertEqual(self.theorem["lambda_minus_delta_residual"], "0")
        self.assertTrue(self.theorem["two_node_symbolic_solution_matches_p2503_delta"])

    def test_finite_pair_recovery(self) -> None:
        self.assertEqual(self.theorem["finite_pair_recovery_row_count"], 55)
        self.assertLess(float(self.theorem["max_abs_lambda_minus_delta"]), 1e-75)
        self.assertLess(float(self.theorem["max_abs_beta0_minus_one"]), 1e-75)
        self.assertLess(float(self.theorem["max_abs_pair_reconstruction_residual"]), 1e-70)
        self.assertTrue(self.theorem["all_pairs_recover_same_delta_within_precision"])
        self.assertTrue(self.theorem["all_pairs_recover_beta0_one_within_precision"])
        self.assertTrue(self.theorem["constant_coefficient_rg_candidate_unique_within_ansatz"])

    def test_negative_controls_and_gatekeepers(self) -> None:
        self.assertTrue(self.theorem["ansatz_uniqueness_not_source_derivation"])
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
        self.assertIn("P2504/S1454", MD.read_text(encoding="utf-8"))
        self.assertIn("P2504/S1454", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2504/S1454", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
