from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2536_s1486_strict_damping_minimal_repair_ideal_certificate import (
    MD,
    OUT,
    build_payload,
    write_markdown,
)


class P2536StrictDampingMinimalRepairIdealTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["strict_damping_minimal_repair_ideal_certificate"]["theorem_export"]
        cls.cert = cls.theorem["strict_damping_minimal_repair_ideal_certificate"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2536")
        self.assertEqual(self.payload["stage_id"], "S1486")
        self.assertIn("MINIMAL_REPAIR_IDEAL", self.payload["status"])
        self.assertIn("ROWWISE_EXACT", self.payload["status"])
        self.assertIn("NO_SOURCE_EXPORT", self.payload["status"])
        self.assertIn("NO_BRIDGE", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_inherited_boundaries_and_repair_counts(self) -> None:
        self.assertTrue(self.theorem["p2535_dual_failure_cover_inherited"])
        self.assertTrue(self.theorem["p2532_strictization_graph_inherited"])
        self.assertEqual(self.theorem["valid_ternary_assignment_count"], 81)
        self.assertEqual(self.theorem["strict_accept_assignment_count"], 1)
        self.assertEqual(self.theorem["failure_assignment_count"], 80)
        self.assertEqual(self.theorem["minimal_repair_set_count"], 80)
        self.assertEqual(self.theorem["minimal_repair_set_size_histogram"], {"1": 8, "2": 24, "3": 32, "4": 16})

    def test_rowwise_minimality_and_subset_counts(self) -> None:
        self.assertTrue(self.theorem["rowwise_minimality_verified"])
        self.assertEqual(self.theorem["proper_repair_subset_count_including_empty"], 544)
        self.assertEqual(self.theorem["proper_repair_subset_count_excluding_empty"], 464)
        self.assertEqual(self.theorem["full_minimal_repair_subset_count"], 80)
        self.assertEqual(self.theorem["all_candidate_repair_subsets_count"], 624)
        for row in self.cert["sample_minimal_repair_rows"]:
            self.assertTrue(row["full_repair_accepts"])
            self.assertEqual(row["accepting_proper_repair_subset_count"], 0)
            self.assertEqual(row["theorem_deficit"], len(row["minimal_repair_actions"]))

    def test_bigrade_and_action_incidence(self) -> None:
        self.assertTrue(self.theorem["repair_bigrade_matches_closed_form"])
        self.assertEqual(len(self.cert["minimal_repair_bigrade_rows"]), 14)
        self.assertEqual(sum(row["row_count"] for row in self.cert["minimal_repair_bigrade_rows"]), 80)
        self.assertEqual(sum(row["total_absent_introduction_actions"] for row in self.cert["minimal_repair_bigrade_rows"]), 108)
        self.assertEqual(sum(row["total_axiom_upgrade_actions"] for row in self.cert["minimal_repair_bigrade_rows"]), 108)
        self.assertEqual(self.theorem["total_absent_source_theorem_introduction_actions"], 108)
        self.assertEqual(self.theorem["total_axiom_to_strict_theorem_upgrade_actions"], 108)
        self.assertEqual(self.theorem["total_minimal_repair_action_incidence"], 216)
        self.assertTrue(self.theorem["uniform_key_repair_action_incidence"])
        self.assertTrue(self.theorem["uniform_absent_and_axiom_action_incidence_per_key"])
        for row in self.cert["key_action_rows"]:
            self.assertEqual(row["absent_source_theorem_introduction_rows"], 27)
            self.assertEqual(row["axiom_to_strict_theorem_upgrade_rows"], 27)
            self.assertEqual(row["total_minimal_repair_action_incidence"], 54)

    def test_negative_controls(self) -> None:
        self.assertTrue(self.theorem["minimal_repair_ideal_certificate_exported"])
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
        self.assertIn("P2536/S1486", MD.read_text(encoding="utf-8"))
        self.assertIn("P2536/S1486", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2536/S1486", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
