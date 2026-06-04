from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2540_s1490_strict_damping_m2_operator_signature_current_premise_obstruction_certificate import (
    MD,
    OUT,
    append_doc_sections,
    build_payload,
    write_markdown,
)


class P2540StrictDampingM2OperatorSignatureCurrentPremiseObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["strict_damping_m2_operator_signature_current_premise_obstruction_certificate"]["theorem_export"]
        cls.cert = cls.theorem["strict_damping_m2_operator_signature_current_premise_obstruction_certificate"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2540")
        self.assertEqual(self.payload["stage_id"], "S1490")
        self.assertIn("M2_OPERATOR_SIGNATURE", self.payload["status"])
        self.assertIn("CURRENT_PREMISE_OBSTRUCTION", self.payload["status"])
        self.assertIn("CURRENT_ROUTE_REFUTED", self.payload["status"])
        self.assertIn("NO_SOURCE_EXPORT", self.payload["status"])
        self.assertIn("NO_BRIDGE", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_precursors_and_extended_finite_audit(self) -> None:
        self.assertTrue(self.theorem["p2509_postulated_functional_wellposedness_inherited"])
        self.assertTrue(self.theorem["p2514_higher_order_nonidentifiability_inherited"])
        self.assertTrue(self.theorem["p2515_m2_signature_identified_but_unsourced_inherited"])
        self.assertTrue(self.theorem["p2516_dual_key_boundary_inherited"])
        self.assertTrue(self.theorem["p2530_four_key_irredundancy_inherited"])
        self.assertTrue(self.theorem["p2539_next_step_recommendation_inherited"])
        self.assertEqual(self.theorem["finite_order_range"], list(range(1, 11)))
        self.assertEqual(self.theorem["finite_basis_degrees"], list(range(7)))
        self.assertTrue(self.theorem["all_orders_have_positive_tangent_gram_on_extended_basis"])

    def test_current_premise_nonentailment_countermodel(self) -> None:
        self.assertEqual(self.theorem["current_source_free_premise_passing_orders"], list(range(1, 11)))
        self.assertEqual(self.theorem["passing_order_count"], 10)
        self.assertEqual(self.theorem["countermodel_pair_orders"], [2, 3])
        self.assertTrue(self.theorem["m2_and_m3_both_pass_current_source_free_premises"])
        self.assertFalse(self.theorem["m2_is_unique_under_current_source_free_premises"])
        self.assertFalse(self.theorem["current_source_free_premises_entail_m2"])
        self.assertTrue(self.theorem["current_premise_nonentailment_of_m2_exported"])
        self.assertTrue(self.theorem["m2_operator_signature_source_route_refuted_for_current_source_free_premises"])

    def test_negative_controls(self) -> None:
        self.assertFalse(self.theorem["multiplicative_character_law_source_exported"])
        self.assertFalse(self.theorem["prime_log_proportionality_source_exported"])
        self.assertFalse(self.theorem["slope_value_or_prime_anchor_source_exported"])
        self.assertFalse(self.theorem["beta_eta_numeric_source_exported"])
        self.assertFalse(self.theorem["m2_operator_signature_source_exported"])
        self.assertFalse(self.theorem["strict_damping_beta_eta_source_exported"])
        self.assertFalse(self.theorem["source_obligation_discharge_exported"])
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
        self.assertIn("P2540/S1490", MD.read_text(encoding="utf-8"))
        self.assertIn("P2540/S1490", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2540/S1490", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
