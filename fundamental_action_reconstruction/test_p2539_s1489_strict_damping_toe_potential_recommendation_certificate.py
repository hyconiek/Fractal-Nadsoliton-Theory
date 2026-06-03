from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2539_s1489_strict_damping_toe_potential_recommendation_certificate import (
    MD,
    NEXT_STEP_RECOMMENDATION,
    OUT,
    build_payload,
    write_markdown,
)


class P2539StrictDampingToEPotentialRecommendationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["strict_damping_toe_potential_recommendation_certificate"]["theorem_export"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2539")
        self.assertEqual(self.payload["stage_id"], "S1489")
        self.assertIn("TOE_POTENTIAL_RECOMMENDATION", self.payload["status"])
        self.assertIn("ZERO_TOE_GATE_DELTA", self.payload["status"])
        self.assertIn("NEXT_SOURCE_THEOREM", self.payload["status"])
        self.assertIn("NO_SOURCE_EXPORT", self.payload["status"])
        self.assertIn("NO_BRIDGE", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_inherited_toe_frontier_and_gate_delta(self) -> None:
        self.assertTrue(self.theorem["p2538_rewrite_normalization_inherited"])
        self.assertTrue(self.theorem["p2421_toe_prime_implicant_inherited"])
        self.assertEqual(len(self.theorem["toe_gate_names"]), 7)
        self.assertEqual(self.theorem["current_toe_readiness_fraction"], "1/128")
        self.assertEqual(self.theorem["current_toe_repair_distance"], 5)
        self.assertEqual(self.theorem["strict_damping_bookkeeping_toe_gate_delta"], 0)
        self.assertFalse(self.theorem["bookkeeping_changes_toe_readiness"])

    def test_local_subkey_table_and_scenarios(self) -> None:
        self.assertEqual(self.theorem["local_subkey_truth_table_row_count"], 16)
        self.assertEqual(self.theorem["local_subkey_full_rows"], 1)
        self.assertEqual(self.theorem["local_subkey_proper_rows"], 15)
        self.assertTrue(self.theorem["local_subkey_no_direct_toe_gate_flip_verified"])
        self.assertEqual(self.theorem["scenario_count"], 4)
        scenarios = {row["scenario"]: row for row in self.theorem["toe_potential_scenarios"]}
        self.assertFalse(scenarios["after_p2536_p2538_bookkeeping_only"]["toe_ready"])
        self.assertEqual(scenarios["after_p2536_p2538_bookkeeping_only"]["repair_distance_to_toe"], 5)
        self.assertFalse(self.theorem["source_gate_alone_closes_toe"])
        self.assertTrue(self.theorem["conditional_all_missing_gates_close_toe"])

    def test_recommendation_and_negative_controls(self) -> None:
        self.assertEqual(self.theorem["recommended_next_honest_step"], NEXT_STEP_RECOMMENDATION)
        self.assertIn("zero_current_toe_gate_delta", self.theorem["toe_potential_assessment"])
        self.assertTrue(self.theorem["toe_potential_recommendation_certificate_exported"])
        self.assertFalse(self.theorem["axiom_promotion_to_strict_exported"])
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
        self.assertIn("P2539/S1489", MD.read_text(encoding="utf-8"))
        self.assertIn("P2539/S1489", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2539/S1489", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
