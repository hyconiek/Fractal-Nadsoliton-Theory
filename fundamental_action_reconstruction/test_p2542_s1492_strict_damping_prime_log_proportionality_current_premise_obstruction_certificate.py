from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2542_s1492_strict_damping_prime_log_proportionality_current_premise_obstruction_certificate import (
    MD,
    OUT,
    append_doc_sections,
    build_payload,
    write_markdown,
)


class P2542StrictDampingPrimeLogProportionalityCurrentPremiseObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["strict_damping_prime_log_proportionality_current_premise_obstruction_certificate"]["theorem_export"]
        cls.cert = cls.theorem["strict_damping_prime_log_proportionality_current_premise_obstruction_certificate"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2542")
        self.assertEqual(self.payload["stage_id"], "S1492")
        self.assertIn("PRIME_LOG_PROPORTIONALITY", self.payload["status"])
        self.assertIn("CURRENT_PREMISE_OBSTRUCTION", self.payload["status"])
        self.assertIn("CURRENT_MONOID_ROUTE_REFUTED", self.payload["status"])
        self.assertIn("NO_SOURCE_EXPORT", self.payload["status"])
        self.assertIn("NO_BRIDGE", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_precursors_inherited(self) -> None:
        self.assertTrue(self.theorem["p2526_prime_character_nullity_inherited"])
        self.assertTrue(self.theorem["p2527_slope_line_subkey_inherited"])
        self.assertTrue(self.theorem["p2530_four_key_irredundancy_inherited"])
        self.assertTrue(self.theorem["p2539_next_step_recommendation_inherited"])
        self.assertTrue(self.theorem["p2541_multiplicative_route_obstruction_inherited"])

    def test_prime_character_countermodels(self) -> None:
        self.assertEqual(self.theorem["sample_prime_character_row_count"], 7)
        self.assertTrue(self.theorem["all_sample_rows_are_unital_multiplicative_characters"])
        self.assertGreaterEqual(self.theorem["countermodel_count"], 4)
        self.assertTrue(self.theorem["current_unital_multiplicative_premises_do_not_entail_prime_log_proportionality"])
        witness = self.theorem["chosen_countermodel_unit_p2_only"]
        self.assertEqual(witness["sample_name"], "unit_p2_only_character")
        self.assertTrue(witness["multiplicative_character_accepts"])
        self.assertFalse(witness["prime_log_proportionality_accepts"])
        self.assertGreater(witness["normalized_ratio_spread"], 0)

    def test_ratio_collapse_and_slope_overgeneration(self) -> None:
        self.assertEqual(self.theorem["exact_ratio_constraint_rank"], 4)
        self.assertEqual(self.theorem["exact_ratio_constraint_nullity"], 1)
        self.assertTrue(self.theorem["strict_log_slope_is_one_proportional_member"])
        self.assertTrue(self.theorem["nonstrict_log_slope_also_proportional"])
        self.assertTrue(self.theorem["prime_log_proportionality_overgenerates_slope_line"])
        proportional_rows = [row for row in self.cert["prime_character_rows"] if row["prime_log_proportionality_accepts"]]
        self.assertGreaterEqual(len(proportional_rows), 2)

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
        self.assertIn("P2542/S1492", MD.read_text(encoding="utf-8"))
        self.assertIn("P2542/S1492", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2542/S1492", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
