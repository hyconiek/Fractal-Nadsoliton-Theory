from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2524_s1474_strict_damping_affine_consistency_continuum_nonidentifiability_certificate import (
    MD,
    OUT,
    build_payload,
    write_markdown,
)


class P2524StrictDampingAffineConsistencyContinuumNonidentifiabilityTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["strict_damping_affine_consistency_continuum_nonidentifiability_certificate"]["theorem_export"]
        cls.cert = cls.theorem["strict_damping_affine_consistency_continuum_nonidentifiability_certificate"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2524")
        self.assertEqual(self.payload["stage_id"], "S1474")
        self.assertIn("AFFINE_CONSISTENCY_CONTINUUM_NONIDENTIFIABILITY", self.payload["status"])
        self.assertIn("NO_SOURCE_EXPORT", self.payload["status"])
        self.assertIn("NO_BRIDGE", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_affine_consistency_continuum(self) -> None:
        self.assertTrue(self.theorem["p2523_basis_independent_consistency_inherited"])
        self.assertEqual(self.theorem["candidate_grid_row_count"], 35)
        self.assertEqual(self.theorem["affine_consistency_accepting_row_count"], 35)
        self.assertEqual(self.theorem["strict_numeric_target_row_count"], 1)
        self.assertTrue(self.theorem["all_candidate_grid_rows_pass_affine_consistency"])
        self.assertTrue(self.theorem["strict_target_is_one_member_of_affine_consistency_continuum"])
        self.assertTrue(all(row["basis_independent_affine_consistency_accepts"] for row in self.cert["affine_consistency_rows"]))

    def test_numeric_filters_are_separate(self) -> None:
        audit = self.cert["endpoint_anchor_filter_audit"]
        self.assertTrue(audit["left_normalization_alone_leaves_slope_continuum_on_grid"])
        self.assertTrue(audit["slope_value_alone_leaves_intercept_continuum_on_grid"])
        self.assertTrue(audit["both_numeric_filters_unique_strict_target"])
        self.assertEqual(audit["left_normalization_filter_count"], 7)
        self.assertEqual(audit["slope_value_filter_count"], 5)
        self.assertEqual(audit["both_numeric_filters_count"], 1)

    def test_negative_controls(self) -> None:
        self.assertTrue(self.theorem["basis_independent_affine_consistency_nonidentifiability_exported"])
        self.assertFalse(self.theorem["node_data_source_exported"])
        self.assertFalse(self.theorem["left_normalization_source_exported"])
        self.assertFalse(self.theorem["slope_value_source_exported"])
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
        self.assertIn("P2524/S1474", MD.read_text(encoding="utf-8"))
        self.assertIn("P2524/S1474", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2524/S1474", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
