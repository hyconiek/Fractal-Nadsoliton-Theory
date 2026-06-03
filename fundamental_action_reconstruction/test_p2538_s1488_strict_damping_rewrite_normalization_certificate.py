from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2538_s1488_strict_damping_rewrite_normalization_certificate import (
    MD,
    OUT,
    build_payload,
    write_markdown,
)


class P2538StrictDampingRewriteNormalizationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["strict_damping_rewrite_normalization_certificate"]["theorem_export"]
        cls.cert = cls.theorem["strict_damping_rewrite_normalization_certificate"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2538")
        self.assertEqual(self.payload["stage_id"], "S1488")
        self.assertIn("REWRITE_NORMALIZATION", self.payload["status"])
        self.assertIn("FINITE_NEWMAN", self.payload["status"])
        self.assertIn("UNIQUE_NORMAL_FORM", self.payload["status"])
        self.assertIn("NO_SOURCE_EXPORT", self.payload["status"])
        self.assertIn("NO_BRIDGE", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_inherited_boundaries_and_rewrite_graph(self) -> None:
        self.assertTrue(self.theorem["p2537_repair_confluence_cube_inherited"])
        self.assertTrue(self.theorem["p2532_strictization_graph_inherited"])
        self.assertEqual(self.theorem["rewrite_vertex_count"], 81)
        self.assertEqual(self.theorem["rewrite_edge_count"], 216)
        expected_hist = {"0": 1, "1": 8, "2": 24, "3": 32, "4": 16}
        self.assertEqual(self.theorem["rank_histogram"], expected_hist)
        self.assertEqual(self.theorem["distance_to_terminal_histogram"], expected_hist)
        self.assertTrue(self.theorem["rank_equals_shortest_distance_to_terminal"])
        self.assertTrue(self.theorem["rank_drop_by_one_on_every_edge"])

    def test_termination_and_unique_terminal(self) -> None:
        self.assertEqual(self.theorem["terminal_vertex_count"], 1)
        self.assertTrue(self.theorem["all_vertices_reach_unique_terminal"])
        self.assertEqual(self.theorem["strongly_connected_component_count"], 81)
        self.assertEqual(self.theorem["nontrivial_strongly_connected_component_count"], 0)
        self.assertTrue(self.theorem["rewrite_graph_acyclic_by_rank"])
        self.assertTrue(self.theorem["global_unique_normal_form_verified"])

    def test_critical_pair_basis_and_finite_newman(self) -> None:
        self.assertEqual(self.theorem["critical_pair_count"], 216)
        self.assertEqual(self.theorem["critical_pair_typed_counts"], {
            "absent_source_theorem_introduction+absent_source_theorem_introduction": 54,
            "absent_source_theorem_introduction+axiom_to_strict_theorem_upgrade": 108,
            "axiom_to_strict_theorem_upgrade+axiom_to_strict_theorem_upgrade": 54,
        })
        self.assertTrue(self.theorem["all_critical_pairs_join_in_one_step_each"])
        self.assertTrue(self.theorem["finite_newman_conditions_verified"])
        for pair in self.cert["sample_critical_pairs"]:
            self.assertTrue(pair["converges_in_one_more_step_each"])
            self.assertEqual(pair["left_then_right_join"], pair["right_then_left_join"])

    def test_negative_controls(self) -> None:
        self.assertTrue(self.theorem["rewrite_normalization_certificate_exported"])
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
        self.assertIn("P2538/S1488", MD.read_text(encoding="utf-8"))
        self.assertIn("P2538/S1488", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2538/S1488", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
