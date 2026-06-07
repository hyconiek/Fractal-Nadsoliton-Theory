from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2546_s1496_strict_damping_identity_action_conditional_propagation_certificate import (
    MD,
    OUT,
    append_doc_sections,
    build_payload,
    write_markdown,
)


class P2546StrictDampingIdentityActionConditionalPropagationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["strict_damping_identity_action_conditional_propagation_certificate"]["theorem_export"]
        cls.cert = cls.theorem["strict_damping_identity_action_conditional_propagation_certificate"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2546")
        self.assertEqual(self.payload["stage_id"], "S1496")
        self.assertIn("IDENTITY_ACTION", self.payload["status"])
        self.assertIn("CONDITIONAL_PROPAGATION", self.payload["status"])
        self.assertIn("NO_IDENTITY_SOURCE_EXPORT", self.payload["status"])
        self.assertIn("NO_BRIDGE", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_precursors_inherited(self) -> None:
        self.assertTrue(self.theorem["p2530_four_key_irredundancy_inherited"])
        self.assertTrue(self.theorem["p2541_unital_to_multiplicative_equivalence_inherited"])
        self.assertTrue(self.theorem["p2542_prime_log_obstruction_inherited"])
        self.assertTrue(self.theorem["p2543_slope_value_obstruction_inherited"])
        self.assertTrue(self.theorem["p2544_four_key_blocker_inherited"])
        self.assertTrue(self.theorem["p2545_unital_obstruction_inherited"])

    def test_exact_identity_action_audit(self) -> None:
        exact = self.theorem["identity_action_exact_linear_audit"]
        self.assertEqual(exact["solution_form"], "b=0 with free slope a")
        self.assertEqual(exact["rank"], 1)
        self.assertEqual(exact["nullity_after_unit_action"], 1)
        self.assertTrue(exact["unit_action_entails_unital_normalization"])
        self.assertTrue(exact["unit_action_does_not_select_slope"])
        self.assertTrue(self.theorem["identity_action_equivalent_to_y1_zero_on_reused_grid"])
        self.assertTrue(self.theorem["identity_action_rejects_all_p2545_y1_countermodels"])

    def test_conditional_propagation_is_one_key_only(self) -> None:
        self.assertEqual(self.theorem["candidate_rows_reused_from_p2545"], 25)
        self.assertEqual(self.theorem["identity_action_accepting_row_count"], 5)
        self.assertEqual(self.theorem["identity_action_rejecting_row_count"], 20)
        self.assertEqual(self.theorem["conditional_missing_source_key_delta"], 1)
        conditional = self.theorem["conditional_assignment_with_strict_identity_action_source"]
        self.assertTrue(conditional["assignment"]["multiplicative_character_law_source"])
        self.assertFalse(conditional["assignment"]["prime_log_proportionality_source"])
        self.assertFalse(conditional["assignment"]["slope_value_or_prime_anchor_source"])
        self.assertFalse(conditional["assignment"]["m2_operator_signature_source"])
        self.assertEqual(conditional["missing_source_key_count"], 3)
        self.assertFalse(conditional["beta_eta_numeric_source_by_assignment"])
        self.assertFalse(conditional["strict_damping_beta_eta_source_by_assignment"])

    def test_negative_controls_and_recommendation(self) -> None:
        self.assertTrue(self.theorem["identity_action_conditional_propagation_exported"])
        self.assertFalse(self.theorem["strict_identity_action_source_exported"])
        self.assertFalse(self.theorem["unital_monoid_normalization_source_exported"])
        self.assertTrue(self.theorem["multiplicative_character_law_source_exported_conditionally_only"])
        self.assertFalse(self.theorem["prime_log_proportionality_source_exported"])
        self.assertFalse(self.theorem["slope_value_or_prime_anchor_source_exported"])
        self.assertFalse(self.theorem["beta_eta_numeric_source_exported"])
        self.assertFalse(self.theorem["m2_operator_signature_source_exported"])
        self.assertFalse(self.theorem["strict_damping_beta_eta_source_exported"])
        self.assertFalse(self.theorem["source_obligation_discharge_exported"])
        self.assertFalse(self.theorem["full_bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["selector_closure_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_certificate"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertIn("derive the unit law from nadsoliton dynamics", self.theorem["honest_next_step_recommendation"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

    def test_rg_audit_and_documentation_markers(self) -> None:
        audit = self.payload["rg_non_duplication_audit"]
        self.assertEqual(audit["tool"], "rg")
        self.assertIn("intended_research_nonduplication", audit["patterns"])
        self.assertGreaterEqual(audit["patterns"]["new_packet"]["count"], 0)
        self.assertIn("P2546/S1496", MD.read_text(encoding="utf-8"))
        self.assertIn("P2546/S1496", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2546/S1496", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
