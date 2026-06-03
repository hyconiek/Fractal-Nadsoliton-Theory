from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2520_s1470_strict_damping_endpoint_anchor_subkey_lattice_certificate import (
    MD,
    OUT,
    build_payload,
    write_markdown,
)


class P2520StrictDampingEndpointAnchorSubkeyLatticeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["strict_damping_endpoint_anchor_subkey_lattice_certificate"]["theorem_export"]
        cls.cert = cls.theorem["strict_damping_endpoint_anchor_subkey_lattice_certificate"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2520")
        self.assertEqual(self.payload["stage_id"], "S1470")
        self.assertIn("ENDPOINT_ANCHOR_SUBKEY_LATTICE", self.payload["status"])
        self.assertIn("SOURCE_OBLIGATION_REFINEMENT", self.payload["status"])
        self.assertIn("NO_SOURCE_EXPORT", self.payload["status"])
        self.assertIn("NO_BRIDGE", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_subkey_rank_lattice(self) -> None:
        self.assertTrue(self.theorem["p2519_endpoint_anchor_acceptance_inherited"])
        self.assertEqual(self.theorem["subkey_lattice_row_count"], 4)
        self.assertEqual(self.theorem["accepted_row_count"], 1)
        self.assertEqual(self.theorem["proper_subkey_obstruction_count"], 3)
        self.assertTrue(self.theorem["only_both_subkeys_accept"])
        self.assertTrue(self.theorem["left_subkey_alone_is_proper_obstruction"])
        self.assertTrue(self.theorem["right_subkey_alone_is_proper_obstruction"])
        rows = self.cert["subkey_rank_lattice"]["rows"]
        accepted = [row for row in rows if row["unique_numeric_beta_eta_target"]]
        self.assertEqual(len(accepted), 1)
        self.assertEqual(set(accepted[0]["active_anchors"]), set(self.cert["subkey_rank_lattice"]["anchors"]))

    def test_candidate_grid_audit(self) -> None:
        audit = self.cert["candidate_pair_audit"]
        self.assertEqual(audit["grid_row_count"], 35)
        self.assertEqual(audit["accepted_by_both_count"], 1)
        self.assertTrue(audit["only_strict_pair_accepted_by_both"])
        accepted = audit["accepted_by_both_unique_pair"]
        self.assertEqual(accepted["intercept_log_beta_candidate"], "0")
        self.assertEqual(accepted["slope_delta_candidate"], "4/5")
        self.assertEqual(accepted["eta_candidate_if_slope_delta"], "9/5")

    def test_negative_controls(self) -> None:
        self.assertTrue(self.theorem["numeric_source_refinement_exported"])
        self.assertFalse(self.theorem["beta_normalization_left_anchor_source_exported"])
        self.assertFalse(self.theorem["right_endpoint_value_anchor_source_exported"])
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
        self.assertIn("P2520/S1470", MD.read_text(encoding="utf-8"))
        self.assertIn("P2520/S1470", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2520/S1470", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
