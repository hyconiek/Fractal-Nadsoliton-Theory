from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2523_s1473_strict_damping_pairwise_secant_consistency_certificate import (
    MD,
    OUT,
    build_payload,
    write_markdown,
)


class P2523StrictDampingPairwiseSecantConsistencyTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["strict_damping_pairwise_secant_consistency_certificate"]["theorem_export"]
        cls.cert = cls.theorem["strict_damping_pairwise_secant_consistency_certificate"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2523")
        self.assertEqual(self.payload["stage_id"], "S1473")
        self.assertIn("PAIRWISE_SECANT_CONSISTENCY", self.payload["status"])
        self.assertIn("BASIS_INDEPENDENT_CHECK", self.payload["status"])
        self.assertIn("NO_SOURCE_EXPORT", self.payload["status"])
        self.assertIn("NO_BRIDGE", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_pairwise_and_triangle_consistency(self) -> None:
        self.assertTrue(self.theorem["p2522_anchor_basis_equivalence_inherited"])
        self.assertEqual(self.theorem["pair_count"], 55)
        self.assertEqual(self.theorem["triangle_count"], 165)
        self.assertTrue(self.theorem["all_pairwise_secants_match_delta"])
        self.assertTrue(self.theorem["all_pairwise_etas_match_9_over_5"])
        self.assertTrue(self.theorem["all_triangle_additive_cocycles_zero"])
        self.assertTrue(self.theorem["all_triangle_secants_match_delta"])
        self.assertTrue(all(row["secant_delta_matches_4_over_5"] for row in self.cert["pairwise_secant_rows"]))
        self.assertTrue(all(row["additive_y_cocycle_zero"] for row in self.cert["triangle_cocycle_rows"]))

    def test_affine_projection(self) -> None:
        design = self.cert["affine_design_projection_audit"]
        self.assertTrue(design["normal_matrix_determinant_positive"])
        self.assertTrue(self.theorem["projected_intercept_matches_zero"])
        self.assertTrue(self.theorem["projected_delta_matches_4_over_5"])
        self.assertTrue(self.theorem["projected_eta_matches_9_over_5"])
        self.assertTrue(self.theorem["projection_residual_zero"])

    def test_negative_controls(self) -> None:
        self.assertTrue(self.theorem["basis_independent_affine_consistency_exported"])
        self.assertFalse(self.theorem["node_data_source_exported"])
        self.assertFalse(self.theorem["anchor_basis_source_exported"])
        self.assertFalse(self.theorem["beta_eta_numeric_source_exported"])
        self.assertFalse(self.theorem["m2_operator_signature_source_exported"])
        self.assertFalse(self.theorem["strict_damping_beta_eta_source_exported"])
        self.assertFalse(self.theorem["damping_compression_bridge_component_ready"])
        self.assertFalse(self.theorem["full_bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["selector_closure_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_certificate"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

    def test_rg_audit_and_documentation_markers(self) -> None:
        audit = self.payload["rg_non_duplication_audit"]
        self.assertEqual(audit["tool"], "rg")
        self.assertIn("new_packet", audit["patterns"])
        self.assertGreaterEqual(audit["patterns"]["new_packet"]["count"], 0)
        self.assertIn("P2523/S1473", MD.read_text(encoding="utf-8"))
        self.assertIn("P2523/S1473", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2523/S1473", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
