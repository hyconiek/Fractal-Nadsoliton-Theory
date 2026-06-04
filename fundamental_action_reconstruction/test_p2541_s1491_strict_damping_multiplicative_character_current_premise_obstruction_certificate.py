from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2541_s1491_strict_damping_multiplicative_character_current_premise_obstruction_certificate import (
    MD,
    OUT,
    append_doc_sections,
    build_payload,
    write_markdown,
)


class P2541StrictDampingMultiplicativeCharacterCurrentPremiseObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["strict_damping_multiplicative_character_current_premise_obstruction_certificate"]["theorem_export"]
        cls.cert = cls.theorem["strict_damping_multiplicative_character_current_premise_obstruction_certificate"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2541")
        self.assertEqual(self.payload["stage_id"], "S1491")
        self.assertIn("MULTIPLICATIVE_CHARACTER", self.payload["status"])
        self.assertIn("CURRENT_PREMISE_OBSTRUCTION", self.payload["status"])
        self.assertIn("CURRENT_AFFINE_ROUTE_REFUTED", self.payload["status"])
        self.assertIn("NO_SOURCE_EXPORT", self.payload["status"])
        self.assertIn("NO_BRIDGE", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_precursors_inherited(self) -> None:
        self.assertTrue(self.theorem["p2524_affine_continuum_inherited"])
        self.assertTrue(self.theorem["p2525_conditional_beta_normalization_inherited"])
        self.assertTrue(self.theorem["p2526_prime_character_nullity_inherited"])
        self.assertTrue(self.theorem["p2530_four_key_irredundancy_inherited"])
        self.assertTrue(self.theorem["p2539_next_step_recommendation_inherited"])
        self.assertTrue(self.theorem["p2540_m2_route_refutation_inherited"])

    def test_affine_route_obstruction_and_unital_equivalence(self) -> None:
        self.assertEqual(self.theorem["candidate_grid_row_count"], 35)
        self.assertEqual(self.theorem["multiplicative_accepting_count"], 7)
        self.assertGreater(self.theorem["affine_countermodel_count"], 0)
        self.assertTrue(self.theorem["current_affine_premises_do_not_entail_multiplicative_law"])
        witness = self.theorem["chosen_affine_countermodel_passing_affine_but_failing_multiplicativity"]
        self.assertEqual(witness["intercept_log_beta"], "1/2")
        self.assertEqual(witness["slope_delta"], "4/5")
        self.assertTrue(witness["affine_consistency_premise_accepts"])
        self.assertFalse(witness["multiplicative_character_accepts"])
        self.assertTrue(self.theorem["all_multiplicative_rows_are_unital"])
        self.assertTrue(self.theorem["all_unital_affine_rows_are_multiplicative"])
        self.assertTrue(self.theorem["multiplicative_law_equivalent_to_unital_left_normalization_inside_affine_family"])

    def test_multiplicative_law_still_overgenerates_slope(self) -> None:
        self.assertTrue(self.theorem["multiplicative_law_overgenerates_slope_family"])
        self.assertEqual(self.theorem["multiplicative_nonstrict_slope_count"], 6)
        slope_witness = self.theorem["chosen_multiplicative_nonstrict_slope_witness"]
        self.assertEqual(slope_witness["intercept_log_beta"], "0")
        self.assertEqual(slope_witness["slope_delta"], "1/2")
        self.assertTrue(slope_witness["multiplicative_character_accepts"])
        self.assertFalse(slope_witness["is_strict_numeric_target"])

    def test_negative_controls(self) -> None:
        self.assertTrue(self.theorem["current_premise_obstruction_exported"])
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
        self.assertIn("P2541/S1491", MD.read_text(encoding="utf-8"))
        self.assertIn("P2541/S1491", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2541/S1491", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
