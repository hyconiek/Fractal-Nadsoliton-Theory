from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2600_s1550_strict_damping_post_m2_residual_source_matrix import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2600StrictDampingPostM2ResidualSourceMatrixTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["strict_damping_post_m2_residual_source_matrix"]["theorem_export"]
        cls.matrix = cls.theorem["strict_damping_post_m2_residual_source_matrix"]

    def test_identity_and_inheritance(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2600")
        self.assertEqual(self.payload["stage_id"], "S1550")
        self.assertIn("POST_M2_RESIDUAL_SOURCE_MATRIX", self.payload["status"])
        self.assertTrue(self.matrix["p2530_four_key_irredundancy_inherited"])
        self.assertTrue(self.matrix["p2547_post_identity_residual_trikey_inherited"])
        self.assertTrue(self.matrix["p2599_m2_source_inherited"])
        self.assertTrue(self.theorem["m2_operator_signature_source_exported"])

    def test_residual_truth_table(self) -> None:
        self.assertEqual(self.matrix["residual_truth_table_row_count"], 8)
        self.assertEqual(self.matrix["residual_accepting_row_count"], 1)
        self.assertTrue(self.matrix["all_single_residual_omissions_reject_strict_damping"])
        self.assertEqual(self.matrix["post_m2_missing_source_key_count_by_current_assignment"], 3)
        accepting = self.matrix["residual_accepting_row"]
        self.assertTrue(accepting["strict_damping_beta_eta_source_accepts"])
        self.assertTrue(all(accepting["assignment"].values()))
        self.assertFalse(self.theorem["strict_damping_beta_eta_source_exported"])
        self.assertTrue(self.theorem["strict_damping_beta_eta_source_remains_blocked_by_current_assignment"])

    def test_scope_guards_and_docs(self) -> None:
        self.assertFalse(self.theorem["multiplicative_character_law_source_exported"])
        self.assertFalse(self.theorem["prime_log_proportionality_source_exported"])
        self.assertFalse(self.theorem["slope_value_or_prime_anchor_source_exported"])
        self.assertFalse(self.theorem["beta_eta_numeric_source_exported"])
        self.assertFalse(self.theorem["bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_certificate"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2600/S1550", MD.read_text(encoding="utf-8"))
        self.assertIn("P2600/S1550", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2600/S1550", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
