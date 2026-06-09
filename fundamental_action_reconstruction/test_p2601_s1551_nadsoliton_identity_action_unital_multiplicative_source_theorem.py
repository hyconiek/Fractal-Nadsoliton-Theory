from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2601_s1551_nadsoliton_identity_action_unital_multiplicative_source_theorem import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2601NadsolitonIdentityActionUnitalMultiplicativeSourceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["nadsoliton_identity_action_unital_multiplicative_source_theorem"]["theorem_export"]
        cls.derivation = cls.theorem["identity_action_derivation"]
        cls.residual = cls.theorem["post_unital_post_m2_residual_matrix"]

    def test_identity_and_inheritance(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2601")
        self.assertEqual(self.payload["stage_id"], "S1551")
        self.assertIn("UNITAL_MULTIPLICATIVE_SOURCE", self.payload["status"])
        self.assertTrue(self.payload["gatekeeper_checks"]["p2541_equivalence_inherited"])
        self.assertTrue(self.payload["gatekeeper_checks"]["p2546_conditional_identity_route_inherited"])
        self.assertTrue(self.payload["gatekeeper_checks"]["p2600_m2_source_inherited"])

    def test_identity_action_closes_multiplicative_key(self) -> None:
        self.assertTrue(self.theorem["hydrodynamic_identity_action_source_exported"])
        self.assertTrue(self.theorem["unital_left_normalization_source_exported"])
        self.assertTrue(self.theorem["multiplicative_character_law_source_exported"])
        self.assertEqual(self.derivation["rg_time_tau_log_lambda_at_identity"], "0")
        self.assertEqual(self.derivation["damping_amplitude_at_identity"], "y_1 = integral_0^0 gamma(tau) d tau = 0")
        self.assertEqual(self.derivation["candidate_row_count"], 9)
        self.assertTrue(self.derivation["all_unital_rows_are_multiplicative"])
        self.assertTrue(self.derivation["all_nonunital_rows_fail_multiplicativity"])

    def test_post_unital_post_m2_residual_matrix(self) -> None:
        self.assertEqual(self.residual["remaining_keys_after_m2_and_unital"], [
            "prime_log_proportionality_source",
            "slope_value_or_prime_anchor_source",
        ])
        self.assertEqual(self.residual["residual_truth_table_row_count"], 4)
        self.assertEqual(self.residual["residual_accepting_row_count"], 1)
        self.assertEqual(self.residual["remaining_missing_source_key_count_after_p2601"], 2)
        self.assertFalse(self.residual["strict_damping_beta_eta_source_accepts_after_p2601_current_assignment"])
        self.assertFalse(self.theorem["strict_damping_beta_eta_source_exported"])

    def test_scope_guards_and_docs(self) -> None:
        self.assertFalse(self.theorem["prime_log_proportionality_source_exported"])
        self.assertFalse(self.theorem["slope_value_or_prime_anchor_source_exported"])
        self.assertFalse(self.theorem["beta_eta_numeric_source_exported"])
        self.assertFalse(self.theorem["bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_theorem"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2601/S1551", MD.read_text(encoding="utf-8"))
        self.assertIn("P2601/S1551", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2601/S1551", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
