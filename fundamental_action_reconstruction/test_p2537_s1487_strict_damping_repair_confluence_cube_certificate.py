from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2537_s1487_strict_damping_repair_confluence_cube_certificate import (
    MD,
    OUT,
    build_payload,
    write_markdown,
)


class P2537StrictDampingRepairConfluenceCubeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["strict_damping_repair_confluence_cube_certificate"]["theorem_export"]
        cls.cert = cls.theorem["strict_damping_repair_confluence_cube_certificate"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2537")
        self.assertEqual(self.payload["stage_id"], "S1487")
        self.assertIn("REPAIR_CONFLUENCE_CUBE", self.payload["status"])
        self.assertIn("UNIQUE_NORMAL_FORM", self.payload["status"])
        self.assertIn("NO_SOURCE_EXPORT", self.payload["status"])
        self.assertIn("NO_BRIDGE", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_inherited_boundaries_and_cube_counts(self) -> None:
        self.assertTrue(self.theorem["p2536_minimal_repair_ideal_inherited"])
        self.assertTrue(self.theorem["p2532_strictization_graph_inherited"])
        self.assertEqual(self.theorem["failure_assignment_count"], 80)
        self.assertEqual(self.theorem["repair_dimension_histogram"], {"1": 8, "2": 24, "3": 32, "4": 16})
        self.assertEqual(self.theorem["repair_cube_vertex_count_total"], 624)
        self.assertEqual(self.theorem["repair_cube_edge_count_total"], 1000)
        self.assertEqual(self.theorem["repair_cube_diamond_square_count_total"], 600)
        self.assertTrue(self.theorem["closed_forms_match_enumeration"])

    def test_confluence_and_shortest_paths(self) -> None:
        self.assertTrue(self.theorem["all_row_cubes_confluent"])
        self.assertEqual(self.theorem["shortest_repair_path_count_total"], 632)
        self.assertEqual(self.theorem["shortest_repair_path_edge_traversal_total"], 2216)
        for row in self.cert["sample_confluence_rows"]:
            self.assertTrue(row["all_shortest_paths_have_same_terminal"])
            self.assertTrue(row["terminal_strict_accept"])
            self.assertTrue(row["all_diamonds_commute"])
            self.assertEqual(row["vertex_count"], 2 ** row["repair_dimension"])

    def test_typed_diamond_counts(self) -> None:
        self.assertEqual(self.theorem["diamond_square_typed_counts"], {
            "absent_source_theorem_introduction+absent_source_theorem_introduction": 150,
            "absent_source_theorem_introduction+axiom_to_strict_theorem_upgrade": 300,
            "axiom_to_strict_theorem_upgrade+axiom_to_strict_theorem_upgrade": 150,
        })
        self.assertEqual(sum(self.theorem["diamond_square_typed_counts"].values()), 600)

    def test_negative_controls(self) -> None:
        self.assertTrue(self.theorem["repair_confluence_cube_certificate_exported"])
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
        self.assertIn("P2537/S1487", MD.read_text(encoding="utf-8"))
        self.assertIn("P2537/S1487", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2537/S1487", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
