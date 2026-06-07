from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2547_s1497_strict_damping_post_identity_residual_trikey_certificate import (
    MD,
    OUT,
    append_doc_sections,
    build_payload,
    write_markdown,
)


class P2547StrictDampingPostIdentityResidualTriKeyTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["strict_damping_post_identity_residual_trikey_certificate"]["theorem_export"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2547")
        self.assertEqual(self.payload["stage_id"], "S1497")
        self.assertIn("POST_IDENTITY_RESIDUAL_TRIKEY", self.payload["status"])
        self.assertIn("NO_RESIDUAL_SOURCE_EXPORT", self.payload["status"])
        self.assertIn("NO_BRIDGE", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_precursors_inherited(self) -> None:
        self.assertTrue(self.theorem["p2530_four_key_irredundancy_inherited"])
        self.assertTrue(self.theorem["p2540_m2_obstruction_inherited"])
        self.assertTrue(self.theorem["p2542_prime_log_obstruction_inherited"])
        self.assertTrue(self.theorem["p2543_slope_obstruction_inherited"])
        self.assertTrue(self.theorem["p2544_four_key_blocker_inherited"])
        self.assertTrue(self.theorem["p2546_identity_conditional_inherited"])

    def test_residual_truth_table(self) -> None:
        self.assertEqual(self.theorem["residual_source_keys"], [
            "prime_log_proportionality_source",
            "slope_value_or_prime_anchor_source",
            "m2_operator_signature_source",
        ])
        self.assertEqual(self.theorem["residual_truth_table_row_count"], 8)
        self.assertEqual(self.theorem["residual_strict_accepting_row_count"], 1)
        accepting = self.theorem["residual_strict_accepting_row"]
        self.assertTrue(all(accepting["residual_assignment"].values()))
        self.assertTrue(accepting["strict_damping_beta_eta_source_accepts"])

    def test_single_omission_irredundancy(self) -> None:
        witnesses = self.theorem["residual_single_omission_witnesses"]
        self.assertEqual(len(witnesses), 3)
        self.assertTrue(self.theorem["all_residual_single_omissions_reject_strict_damping"])
        self.assertEqual(
            sorted(witness["omitted_key"] for witness in witnesses),
            sorted(self.theorem["residual_source_keys"]),
        )
        self.assertTrue(self.theorem["conditional_identity_reduces_missing_count_from_4_to_3"])
        self.assertTrue(self.theorem["identity_source_alone_cannot_export_beta_eta_numeric_source"])
        self.assertTrue(self.theorem["identity_source_alone_cannot_export_strict_damping_beta_eta_source"])

    def test_negative_controls_and_recommendation(self) -> None:
        self.assertTrue(self.theorem["post_identity_residual_trikey_irredundancy_exported"])
        self.assertFalse(self.theorem["prime_log_proportionality_source_exported"])
        self.assertFalse(self.theorem["slope_value_or_prime_anchor_source_exported"])
        self.assertFalse(self.theorem["m2_operator_signature_source_exported"])
        self.assertFalse(self.theorem["beta_eta_numeric_source_exported"])
        self.assertFalse(self.theorem["strict_damping_beta_eta_source_exported"])
        self.assertFalse(self.theorem["source_obligation_discharge_exported"])
        self.assertFalse(self.theorem["full_bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["selector_closure_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_certificate"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertIn("m=2 operator-order selection", self.theorem["recommended_next_honest_step"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

    def test_rg_audit_and_documentation_markers(self) -> None:
        audit = self.payload["rg_non_duplication_audit"]
        self.assertEqual(audit["tool"], "rg")
        self.assertIn("intended_research_nonduplication", audit["patterns"])
        self.assertGreaterEqual(audit["patterns"]["new_packet"]["count"], 0)
        self.assertIn("P2547/S1497", MD.read_text(encoding="utf-8"))
        self.assertIn("P2547/S1497", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2547/S1497", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
