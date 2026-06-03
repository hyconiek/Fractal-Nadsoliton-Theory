from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2519_s1469_strict_damping_biharmonic_endpoint_anchor_acceptance_certificate import (
    MD,
    OUT,
    build_payload,
    write_markdown,
)


class P2519StrictDampingBiharmonicEndpointAnchorAcceptanceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["strict_damping_biharmonic_endpoint_anchor_acceptance_certificate"]["theorem_export"]
        cls.cert = cls.theorem["strict_damping_biharmonic_endpoint_anchor_acceptance_certificate"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2519")
        self.assertEqual(self.payload["stage_id"], "S1469")
        self.assertIn("ENDPOINT_ANCHOR_ACCEPTANCE", self.payload["status"])
        self.assertIn("CONDITIONAL_TARGET_RECOVERY", self.payload["status"])
        self.assertIn("NO_SOURCE_EXPORT", self.payload["status"])
        self.assertIn("NO_BRIDGE", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_endpoint_anchor_solution(self) -> None:
        self.assertTrue(self.theorem["p2518_nonidentifiability_inherited"])
        self.assertEqual(self.theorem["anchor_determinant"], "log(11)")
        self.assertTrue(self.theorem["anchor_determinant_positive"])
        self.assertEqual(self.theorem["solved_slope_delta"], "4/5")
        self.assertEqual(self.theorem["solved_beta"], "1")
        self.assertEqual(self.theorem["solved_eta"], "9/5")
        self.assertTrue(self.theorem["unique_candidate_slope_accepted_by_endpoint_anchors"])
        self.assertTrue(self.theorem["endpoint_anchors_pin_affine_slope_if_anchors_are_sourced"])

    def test_finite_design_and_node_residuals(self) -> None:
        design = self.cert["finite_design_matrix_audit"]
        self.assertTrue(design["normal_matrix_determinant_positive"])
        self.assertTrue(design["least_squares_slope_matches_delta_4_over_5"])
        self.assertTrue(design["all_node_residuals_zero_within_float_tolerance"])
        self.assertTrue(self.theorem["finite_design_all_node_residuals_zero"])
        self.assertTrue(all(row["residual"] == 0.0 for row in self.cert["all_node_reconstruction_rows"]))

    def test_negative_controls(self) -> None:
        self.assertTrue(self.theorem["conditional_beta_eta_target_recovered_from_anchors"])
        self.assertFalse(self.theorem["endpoint_anchor_source_exported"])
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
        self.assertIn("P2519/S1469", MD.read_text(encoding="utf-8"))
        self.assertIn("P2519/S1469", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2519/S1469", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
