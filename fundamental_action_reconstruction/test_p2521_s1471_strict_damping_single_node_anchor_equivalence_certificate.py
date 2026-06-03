from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2521_s1471_strict_damping_single_node_anchor_equivalence_certificate import (
    MD,
    OUT,
    build_payload,
    write_markdown,
)


class P2521StrictDampingSingleNodeAnchorEquivalenceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["strict_damping_single_node_anchor_equivalence_certificate"]["theorem_export"]
        cls.cert = cls.theorem["strict_damping_single_node_anchor_equivalence_certificate"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2521")
        self.assertEqual(self.payload["stage_id"], "S1471")
        self.assertIn("SINGLE_NODE_ANCHOR_EQUIVALENCE", self.payload["status"])
        self.assertIn("PLACEMENT_NONIDENTIFIABILITY", self.payload["status"])
        self.assertIn("NO_SOURCE_EXPORT", self.payload["status"])
        self.assertIn("NO_BRIDGE", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_single_node_anchor_equivalence(self) -> None:
        self.assertTrue(self.theorem["p2520_subkey_lattice_inherited"])
        self.assertEqual(self.theorem["nonzero_anchor_count"], 10)
        self.assertTrue(self.theorem["every_nonzero_node_with_left_anchor_pins_same_delta"])
        self.assertTrue(self.theorem["every_nonzero_node_has_positive_determinant"])
        self.assertTrue(self.theorem["every_anchor_uniquely_accepts_strict_candidate_on_grid"])
        self.assertTrue(self.theorem["anchor_placement_not_identified_by_affine_algebra"])
        self.assertTrue(self.theorem["right_endpoint_d11_is_sufficient_not_unique"])
        self.assertTrue(self.theorem["single_nonzero_node_anchor_equivalence_exported"])

    def test_equivalence_rows(self) -> None:
        lattice = self.cert["single_node_anchor_lattice"]
        self.assertEqual(len(lattice["single_node_equivalence_rows"]), 10)
        self.assertTrue(all(row["determinant_positive"] for row in lattice["single_node_equivalence_rows"]))
        self.assertTrue(all(row["solved_delta_matches_4_over_5"] for row in lattice["single_node_equivalence_rows"]))
        self.assertTrue(all(row["solved_eta_matches_9_over_5"] for row in lattice["single_node_equivalence_rows"]))
        self.assertTrue(all(row["all_nodes_reconstructed_after_pinning"] for row in lattice["single_node_equivalence_rows"]))

    def test_negative_controls(self) -> None:
        self.assertFalse(self.theorem["anchor_placement_source_exported"])
        self.assertFalse(self.theorem["nonzero_node_value_source_exported"])
        self.assertFalse(self.theorem["beta_normalization_left_anchor_source_exported"])
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
        self.assertIn("P2521/S1471", MD.read_text(encoding="utf-8"))
        self.assertIn("P2521/S1471", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2521/S1471", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
