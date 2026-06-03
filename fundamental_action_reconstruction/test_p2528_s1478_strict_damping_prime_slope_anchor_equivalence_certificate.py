from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2528_s1478_strict_damping_prime_slope_anchor_equivalence_certificate import (
    MD,
    OUT,
    build_payload,
    write_markdown,
)


class P2528StrictDampingPrimeSlopeAnchorEquivalenceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["strict_damping_prime_slope_anchor_equivalence_certificate"]["theorem_export"]
        cls.cert = cls.theorem["strict_damping_prime_slope_anchor_equivalence_certificate"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2528")
        self.assertEqual(self.payload["stage_id"], "S1478")
        self.assertIn("PRIME_SLOPE_ANCHOR_EQUIVALENCE", self.payload["status"])
        self.assertIn("CONDITIONAL_ANCHOR_ONLY", self.payload["status"])
        self.assertIn("NO_ANCHOR_SOURCE", self.payload["status"])
        self.assertIn("NO_SLOPE_SOURCE", self.payload["status"])
        self.assertIn("NO_OPERATOR_SOURCE", self.payload["status"])
        self.assertIn("NO_BRIDGE", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_single_prime_anchor_equivalence(self) -> None:
        self.assertTrue(self.theorem["p2527_slope_line_inherited"])
        self.assertEqual(self.theorem["prime_anchor_count"], 5)
        self.assertEqual(self.theorem["slope_candidate_count_per_anchor"], 7)
        self.assertTrue(self.theorem["all_prime_anchors_have_positive_log_determinant"])
        self.assertTrue(self.theorem["all_single_prime_anchors_recover_strict_delta"])
        self.assertTrue(self.theorem["all_single_prime_anchors_recover_strict_eta"])
        self.assertTrue(self.theorem["all_single_prime_anchors_reconstruct_strict_node_line"])
        self.assertTrue(self.theorem["single_prime_anchor_placement_nonunique"])
        self.assertTrue(all(row["recovers_strict_delta"] for row in self.cert["prime_anchor_rows"]))

    def test_candidate_grid_accepts_only_strict_slope(self) -> None:
        grid = self.cert["candidate_grid_audit_rows"]
        self.assertTrue(self.theorem["all_candidate_grid_audits_accept_only_strict_slope"])
        self.assertEqual(len(grid), 5)
        self.assertTrue(all(row["accepted_slope_candidates"] == ["4/5"] for row in grid))

    def test_negative_controls(self) -> None:
        self.assertTrue(self.theorem["conditional_prime_slope_anchor_equivalence_exported"])
        self.assertFalse(self.theorem["prime_anchor_placement_source_exported"])
        self.assertFalse(self.theorem["prime_anchor_value_source_exported"])
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
        self.assertIn("intended_research_nonduplication", audit["patterns"])
        self.assertGreaterEqual(audit["patterns"]["new_packet"]["count"], 0)
        self.assertIn("P2528/S1478", MD.read_text(encoding="utf-8"))
        self.assertIn("P2528/S1478", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2528/S1478", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
