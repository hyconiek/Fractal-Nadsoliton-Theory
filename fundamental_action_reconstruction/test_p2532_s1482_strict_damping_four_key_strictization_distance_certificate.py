from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2532_s1482_strict_damping_four_key_strictization_distance_certificate import (
    MD,
    OUT,
    build_payload,
    write_markdown,
)


class P2532StrictDampingFourKeyStrictizationDistanceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["strict_damping_four_key_strictization_distance_certificate"]["theorem_export"]
        cls.cert = cls.theorem["strict_damping_four_key_strictization_distance_certificate"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2532")
        self.assertEqual(self.payload["stage_id"], "S1482")
        self.assertIn("FOUR_KEY_STRICTIZATION_DISTANCE", self.payload["status"])
        self.assertIn("DEFICIT_STRATIFICATION_AND_GRAPH_ONLY", self.payload["status"])
        self.assertIn("NO_SOURCE_EXPORT", self.payload["status"])
        self.assertIn("NO_BRIDGE", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_strictization_class_counts(self) -> None:
        self.assertTrue(self.theorem["p2531_axiom_boundary_inherited"])
        self.assertEqual(self.theorem["ternary_status_row_count"], 81)
        self.assertEqual(self.theorem["strict_accept_row_count"], 1)
        self.assertEqual(self.theorem["present_non_strict_axiom_augmented_row_count"], 15)
        self.assertEqual(self.theorem["blocked_missing_key_row_count"], 65)
        self.assertEqual(self.cert["class_counts"], {
            "blocked_missing_key": 65,
            "present_non_strict_axiom_augmented": 15,
            "strict_accept": 1,
        })

    def test_theorem_deficit_and_distance_strata(self) -> None:
        self.assertEqual(self.theorem["theorem_deficit_counts"], {
            "0": 1,
            "1": 8,
            "2": 24,
            "3": 32,
            "4": 16,
        })
        self.assertTrue(self.theorem["theorem_deficit_counts_match_closed_form"])
        self.assertEqual(self.theorem["present_axiom_upgrade_distance_counts"], {
            "1": 4,
            "2": 6,
            "3": 4,
            "4": 1,
        })
        self.assertTrue(self.theorem["present_axiom_upgrade_distance_counts_match_binomial"])
        self.assertEqual(self.theorem["blocked_absent_source_gap_counts"], {
            "1": 32,
            "2": 24,
            "3": 8,
            "4": 1,
        })
        self.assertTrue(self.theorem["blocked_absent_source_gap_counts_match_enumeration"])

    def test_nearest_non_strict_and_positive_deficit(self) -> None:
        self.assertTrue(self.theorem["minimal_non_strict_present_rows_are_one_axiom_rows"])
        self.assertEqual(len(self.cert["minimal_present_axiom_upgrade_rows"]), 4)
        self.assertTrue(all(row["axiom_count"] == 1 for row in self.cert["minimal_present_axiom_upgrade_rows"]))
        self.assertEqual(self.theorem["maximal_absent_blocker_row_count"], 1)
        self.assertTrue(self.theorem["all_non_accepting_rows_have_positive_theorem_deficit"])

    def test_one_step_strictization_graph(self) -> None:
        self.assertEqual(self.theorem["one_step_strictization_edge_count"], 216)
        self.assertTrue(self.theorem["one_step_strictization_edge_count_matches_deficit_sum"])
        self.assertEqual(self.theorem["one_step_edge_type_counts"], {
            "absent_source_theorem_introduction": 108,
            "axiom_to_strict_theorem_upgrade": 108,
        })
        self.assertTrue(self.theorem["one_step_edge_type_counts_match_closed_form"])
        self.assertEqual(self.theorem["one_step_edge_source_class_counts"], {
            "blocked_missing_key": 184,
            "present_non_strict_axiom_augmented": 32,
        })
        self.assertTrue(self.theorem["one_step_edge_source_class_counts_match_strata"])
        self.assertEqual(set(self.theorem["upgraded_key_counts"].values()), {54})
        self.assertTrue(self.theorem["upgraded_key_counts_uniform"])
        self.assertEqual(self.theorem["shortest_strictization_path_count_by_deficit"], {
            "0": 1,
            "1": 8,
            "2": 48,
            "3": 192,
            "4": 384,
        })
        self.assertTrue(self.theorem["shortest_strictization_path_counts_match_factorial_closed_form"])
        graph = self.cert["one_step_strictization_graph_summary"]
        self.assertEqual(graph["total_shortest_strictization_path_count"], 633)
        self.assertEqual(len(graph["edge_samples"]), 16)

    def test_negative_controls(self) -> None:
        self.assertTrue(self.theorem["strictization_distance_certificate_exported"])
        self.assertFalse(self.theorem["axiom_promotion_to_strict_exported"])
        self.assertFalse(self.theorem["multiplicative_character_law_source_exported"])
        self.assertFalse(self.theorem["prime_log_proportionality_source_exported"])
        self.assertFalse(self.theorem["slope_value_or_prime_anchor_source_exported"])
        self.assertFalse(self.theorem["beta_eta_numeric_source_exported"])
        self.assertFalse(self.theorem["m2_operator_signature_source_exported"])
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
        self.assertIn("intended_research_nonduplication", audit["patterns"])
        self.assertGreaterEqual(audit["patterns"]["new_packet"]["count"], 0)
        self.assertIn("P2532/S1482", MD.read_text(encoding="utf-8"))
        self.assertIn("P2532/S1482", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2532/S1482", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
