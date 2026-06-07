from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2545_s1495_strict_damping_unital_normalization_current_premise_obstruction_certificate import (
    MD,
    OUT,
    append_doc_sections,
    build_payload,
    write_markdown,
)


class P2545StrictDampingUnitalNormalizationCurrentPremiseObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["strict_damping_unital_normalization_current_premise_obstruction_certificate"]["theorem_export"]
        cls.cert = cls.theorem["strict_damping_unital_normalization_current_premise_obstruction_certificate"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2545")
        self.assertEqual(self.payload["stage_id"], "S1495")
        self.assertIn("UNITAL_NORMALIZATION", self.payload["status"])
        self.assertIn("CURRENT_PREMISE_OBSTRUCTION", self.payload["status"])
        self.assertIn("NO_SOURCE_EXPORT", self.payload["status"])
        self.assertIn("NO_BRIDGE", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_precursors_inherited(self) -> None:
        self.assertTrue(self.theorem["p2521_left_normalization_unsourced_inherited"])
        self.assertTrue(self.theorem["p2525_conditional_multiplicative_to_beta_normalization_inherited"])
        self.assertTrue(self.theorem["p2530_four_key_irredundancy_inherited"])
        self.assertTrue(self.theorem["p2541_unital_equivalence_inherited"])
        self.assertTrue(self.theorem["p2544_no_false_source_blocker_inherited"])

    def test_unital_normalization_countermodels(self) -> None:
        self.assertEqual(self.theorem["candidate_grid_row_count"], 25)
        self.assertEqual(self.theorem["affine_countermodel_count_to_y1_zero"], 20)
        self.assertEqual(self.theorem["strict_slope_countermodel_count_to_y1_zero"], 4)
        self.assertTrue(self.theorem["identity_node_existence_does_not_entail_y1_zero"])
        self.assertTrue(self.theorem["affine_consistency_does_not_entail_y1_zero"])
        self.assertTrue(self.theorem["even_strict_slope_value_does_not_entail_y1_zero"])
        witness = self.theorem["chosen_current_premise_countermodel_to_unital_normalization"]
        self.assertEqual(witness["intercept_log_beta"], "1/2")
        self.assertEqual(witness["slope_delta"], "4/5")
        self.assertTrue(witness["current_affine_consistency_premise_accepts"])
        self.assertTrue(witness["current_domain_contains_unit_node"])
        self.assertFalse(witness["unital_y1_zero_accepts"])
        self.assertFalse(witness["unit_product_law_accepts"])

    def test_equivalence_to_unit_and_product_laws(self) -> None:
        self.assertTrue(self.theorem["unit_product_law_equivalent_to_y1_zero_on_affine_grid"])
        self.assertTrue(self.theorem["full_product_multiplicativity_equivalent_to_y1_zero_on_affine_grid"])
        self.assertTrue(self.theorem["all_product_defects_equal_minus_y1"])
        self.assertTrue(self.theorem["all_unit_defects_equal_minus_y1"])
        self.assertTrue(self.theorem["unital_normalization_is_minimal_missing_source_for_multiplicative_key_inside_affine_family"])

    def test_negative_controls_and_recommendation(self) -> None:
        self.assertTrue(self.theorem["current_premise_obstruction_exported"])
        self.assertTrue(self.theorem["unital_monoid_normalization_current_premise_obstruction_exported"])
        self.assertFalse(self.theorem["unital_monoid_normalization_source_exported"])
        self.assertFalse(self.theorem["multiplicative_character_law_source_exported"])
        self.assertFalse(self.theorem["beta_eta_numeric_source_exported"])
        self.assertFalse(self.theorem["strict_damping_beta_eta_source_exported"])
        self.assertFalse(self.theorem["source_obligation_discharge_exported"])
        self.assertFalse(self.theorem["damping_compression_bridge_component_ready"])
        self.assertFalse(self.theorem["full_bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["selector_closure_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_certificate"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertIn("unit-node normalization theorem", self.theorem["honest_next_step_recommendation"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

    def test_rg_audit_and_documentation_markers(self) -> None:
        audit = self.payload["rg_non_duplication_audit"]
        self.assertEqual(audit["tool"], "rg")
        self.assertIn("intended_research_nonduplication", audit["patterns"])
        self.assertGreaterEqual(audit["patterns"]["new_packet"]["count"], 0)
        self.assertIn("P2545/S1495", MD.read_text(encoding="utf-8"))
        self.assertIn("P2545/S1495", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2545/S1495", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
