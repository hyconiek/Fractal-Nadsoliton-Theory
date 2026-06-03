from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2502_s1452_strict_completion_bridge_minimal_triple_frontier_certificate import (
    MD,
    OUT,
    build_payload,
    write_markdown,
)


class P2502StrictCompletionBridgeMinimalTripleFrontierCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["strict_completion_bridge_minimal_triple_frontier_certificate"]["theorem_export"]
        cls.cert = cls.theorem["bridge_minimal_triple_frontier_certificate"]
        cls.summary = cls.cert["triple_frontier_summary"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2502")
        self.assertEqual(self.payload["stage_id"], "S1452")
        self.assertIn("BRIDGE_MINIMAL_TRIPLE_FRONTIER", self.payload["status"])
        self.assertIn("NO_SOURCE_EXPORT", self.payload["status"])
        self.assertIn("NO_BRIDGE_THEOREM", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_triple_enumeration_and_unique_bridge_triple(self) -> None:
        self.assertEqual(self.theorem["open_atom_count"], 7)
        self.assertEqual(self.theorem["triple_extension_count"], 35)
        self.assertEqual(self.theorem["bridge_closing_triple_count"], 1)
        expected = [
            "strict_damping_beta_eta_source",
            "strict_dynamical_source_for_A_P_D",
            "strict_phase_frequency_source",
        ]
        self.assertEqual(self.theorem["exact_strict_source_triple"], expected)
        self.assertEqual(self.theorem["bridge_closing_triples"][0]["true_atoms"], expected)
        self.assertEqual(self.theorem["bridge_closing_triples"][0]["target_signature_bridge_role_selector_toe"], "1000")
        self.assertTrue(self.theorem["strict_source_triple_closes_only_bridge"])

    def test_role_selector_toe_counts(self) -> None:
        self.assertEqual(self.theorem["role_transfer_closing_triple_count"], 0)
        self.assertEqual(self.theorem["selector_closing_triple_count"], 15)
        self.assertEqual(self.theorem["toe_closing_triple_count"], 0)
        self.assertTrue(self.theorem["role_transfer_still_needs_four_atom_package"])
        self.assertTrue(self.theorem["toe_still_needs_all_seven_atoms"])
        self.assertTrue(self.summary["all_selector_triples_contain_chi11"])

    def test_p2501_curvature_chain_not_promoted_to_source(self) -> None:
        self.assertTrue(self.theorem["p2501_curvature_enclosure_inherited"])
        self.assertTrue(self.theorem["curvature_chain_does_not_discharge_strict_source_atoms"])
        self.assertTrue(self.cert["local_compression_evidence_is_not_an_open_theorem_atom"])
        self.assertFalse(self.cert["adding_local_curvature_chain_alone_changes_frontier_signature"])
        self.assertTrue(self.theorem["strict_source_triple_is_next_bridge_theorem_frontier"])

    def test_negative_controls_and_gatekeepers(self) -> None:
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
        self.assertIn("P2502/S1452", MD.read_text(encoding="utf-8"))
        self.assertIn("P2502/S1452", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2502/S1452", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
