from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2522_s1472_strict_damping_two_node_anchor_basis_equivalence_certificate import (
    MD,
    OUT,
    build_payload,
    write_markdown,
)


class P2522StrictDampingTwoNodeAnchorBasisEquivalenceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["strict_damping_two_node_anchor_basis_equivalence_certificate"]["theorem_export"]
        cls.cert = cls.theorem["strict_damping_two_node_anchor_basis_equivalence_certificate"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2522")
        self.assertEqual(self.payload["stage_id"], "S1472")
        self.assertIn("TWO_NODE_ANCHOR_BASIS_EQUIVALENCE", self.payload["status"])
        self.assertIn("BASIS_NONUNIQUENESS", self.payload["status"])
        self.assertIn("NO_SOURCE_EXPORT", self.payload["status"])
        self.assertIn("NO_BRIDGE", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_all_node_pairs_pin_target(self) -> None:
        self.assertTrue(self.theorem["p2521_single_node_equivalence_inherited"])
        self.assertEqual(self.theorem["node_pair_count"], 55)
        self.assertTrue(self.theorem["every_distinct_node_pair_has_positive_determinant_in_ordered_domain"])
        self.assertTrue(self.theorem["every_distinct_node_pair_derives_beta_normalization"])
        self.assertTrue(self.theorem["every_distinct_node_pair_pins_delta_eta"])
        self.assertTrue(self.theorem["every_distinct_node_pair_reconstructs_all_nodes"])
        rows = self.cert["two_node_basis_rows"]
        self.assertEqual(len(rows), 55)
        self.assertTrue(all(row["solved_intercept_matches_zero"] for row in rows))
        self.assertTrue(all(row["solved_delta_matches_4_over_5"] for row in rows))
        self.assertTrue(all(row["solved_eta_matches_9_over_5"] for row in rows))

    def test_representative_candidate_grid(self) -> None:
        grid = self.cert["representative_candidate_grid"]
        self.assertEqual(grid["row_count"], 140)
        self.assertTrue(grid["every_representative_pair_uniquely_accepts_strict_pair"])
        for row in grid["by_pair_summary"]:
            self.assertEqual(row["accepted_count"], 1)
            self.assertTrue(row["unique_strict_pair_accepted"])
            self.assertEqual(row["accepted_pair"]["intercept_log_beta_candidate"], "0")
            self.assertEqual(row["accepted_pair"]["slope_delta_candidate"], "4/5")
            self.assertEqual(row["accepted_pair"]["eta_candidate_if_slope_delta"], "9/5")

    def test_negative_controls(self) -> None:
        self.assertTrue(self.theorem["anchor_basis_equivalence_exported"])
        self.assertFalse(self.theorem["anchor_basis_source_exported"])
        self.assertFalse(self.theorem["node_pair_value_source_exported"])
        self.assertFalse(self.theorem["anchor_placement_source_exported"])
        self.assertFalse(self.theorem["beta_normalization_left_anchor_source_exported"])
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
        self.assertIn("P2522/S1472", MD.read_text(encoding="utf-8"))
        self.assertIn("P2522/S1472", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2522/S1472", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
