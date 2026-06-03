from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2533_s1483_strict_damping_ternary_generating_polynomial_certificate import (
    MD,
    OUT,
    build_payload,
    coefficient_row,
    write_markdown,
)


class P2533StrictDampingTernaryGeneratingPolynomialTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["strict_damping_ternary_generating_polynomial_certificate"]["theorem_export"]
        cls.cert = cls.theorem["strict_damping_ternary_generating_polynomial_certificate"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2533")
        self.assertEqual(self.payload["stage_id"], "S1483")
        self.assertIn("TERNARY_GENERATING_POLYNOMIAL", self.payload["status"])
        self.assertIn("SYMBOLIC_SOURCE_BOUNDARY_ONLY", self.payload["status"])
        self.assertIn("NO_SOURCE_EXPORT", self.payload["status"])
        self.assertIn("NO_BRIDGE", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_polynomial_normal_forms_and_support(self) -> None:
        self.assertTrue(self.theorem["p2532_strictization_distance_inherited"])
        self.assertEqual(self.theorem["compressed_universe_polynomial"], "(u_absent + a_axiom + s_strict)^4")
        self.assertEqual(self.theorem["strict_accept_polynomial"], "s_strict^4")
        self.assertEqual(self.theorem["present_axiom_augmented_polynomial"], "(a_axiom + s_strict)^4 - s_strict^4")
        self.assertEqual(self.theorem["blocked_missing_key_polynomial"], "(u_absent + a_axiom + s_strict)^4 - (a_axiom + s_strict)^4")
        self.assertEqual(len(self.cert["full_coefficient_rows"]), 15)
        self.assertEqual(self.cert["strict_accept_coefficient_rows"], [coefficient_row(0, 0, 4, 1)])
        self.assertTrue(self.theorem["strict_accept_support_is_single_s4_monomial"])
        self.assertTrue(self.theorem["present_axiom_and_blocked_support_disjoint"])

    def test_polynomial_counts_match_p2532(self) -> None:
        self.assertEqual(self.theorem["class_counts_from_polynomial"], {
            "strict_accept": 1,
            "present_non_strict_axiom_augmented": 15,
            "blocked_missing_key": 65,
        })
        self.assertTrue(self.theorem["class_counts_match_p2532"])
        self.assertEqual(self.theorem["theorem_deficit_counts_from_polynomial"], {
            "0": 1,
            "1": 8,
            "2": 24,
            "3": 32,
            "4": 16,
        })
        self.assertTrue(self.theorem["theorem_deficit_counts_match_p2532"])
        self.assertEqual(self.theorem["present_axiom_upgrade_distance_counts_from_polynomial"], {
            "1": 4,
            "2": 6,
            "3": 4,
            "4": 1,
        })
        self.assertTrue(self.theorem["present_axiom_upgrade_distance_counts_match_p2532"])
        self.assertEqual(self.theorem["blocked_absent_source_gap_counts_from_polynomial"], {
            "1": 32,
            "2": 24,
            "3": 8,
            "4": 1,
        })
        self.assertTrue(self.theorem["blocked_absent_source_gap_counts_match_p2532"])

    def test_derivative_edge_counts(self) -> None:
        self.assertEqual(self.theorem["derivative_absent_axiom_edge_counts_from_polynomial"], {
            "absent_source_theorem_introduction": 108,
            "axiom_to_strict_theorem_upgrade": 108,
        })
        self.assertTrue(self.theorem["derivative_edge_counts_match_p2532"])
        self.assertEqual(set(self.theorem["uniform_key_upgrade_counts_from_polynomial"].values()), {54})
        self.assertTrue(self.theorem["uniform_key_upgrade_counts_match_p2532"])
        edge_counts = self.cert["derivative_edge_counts_from_polynomial"]
        self.assertEqual(edge_counts["total_one_step_edges_from_degree_derivative"], 216)

    def test_negative_controls(self) -> None:
        self.assertTrue(self.theorem["polynomial_certificate_exported"])
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
        self.assertIn("P2533/S1483", MD.read_text(encoding="utf-8"))
        self.assertIn("P2533/S1483", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2533/S1483", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
