from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2516_s1466_strict_damping_dual_key_source_acceptance_matrix import (
    MD,
    OUT,
    build_payload,
    write_markdown,
)


class P2516StrictDampingDualKeySourceAcceptanceMatrixTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["strict_damping_dual_key_source_acceptance_matrix"]["theorem_export"]
        cls.cert = cls.theorem["strict_damping_dual_key_source_acceptance_matrix"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2516")
        self.assertEqual(self.payload["stage_id"], "S1466")
        self.assertIn("DUAL_KEY_SOURCE_ACCEPTANCE", self.payload["status"])
        self.assertIn("NO_SOURCE_EXPORT", self.payload["status"])
        self.assertIn("NO_BRIDGE", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_dual_key_truth_table(self) -> None:
        self.assertEqual(self.theorem["frontier_atom_under_attack"], "strict_damping_beta_eta_source")
        self.assertTrue(self.theorem["p2414_numeric_beta_eta_identified_but_unsourced"])
        self.assertTrue(self.theorem["p2515_m2_operator_signature_identified_but_unsourced"])
        self.assertEqual(
            self.theorem["boolean_normal_form"],
            "strict_damping_beta_eta_source = beta_eta_numeric_source AND m2_operator_signature_source",
        )
        self.assertEqual(self.theorem["accepted_row_count"], 1)
        self.assertEqual(self.theorem["proper_subset_count"], 3)
        self.assertTrue(self.theorem["all_proper_subsets_rejected"])
        self.assertEqual(
            self.theorem["unique_minimal_accepting_set"],
            ["beta_eta_numeric_source", "m2_operator_signature_source"],
        )
        self.assertTrue(self.theorem["numeric_key_alone_insufficient"])
        self.assertTrue(self.theorem["operator_signature_key_alone_insufficient"])
        self.assertTrue(self.theorem["both_keys_required_for_future_source_theorem"])

    def test_negative_controls(self) -> None:
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
        self.assertIn("P2516/S1466", MD.read_text(encoding="utf-8"))
        self.assertIn("P2516/S1466", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2516/S1466", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
