from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2535_s1485_strict_damping_dual_failure_cover_certificate import (
    MD,
    OUT,
    build_payload,
    write_markdown,
)


class P2535StrictDampingDualFailureCoverTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["strict_damping_dual_failure_cover_certificate"]["theorem_export"]
        cls.cert = cls.theorem["strict_damping_dual_failure_cover_certificate"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2535")
        self.assertEqual(self.payload["stage_id"], "S1485")
        self.assertIn("DUAL_FAILURE_COVER", self.payload["status"])
        self.assertIn("ALL_NONSTRICT_KEYS_BLOCK", self.payload["status"])
        self.assertIn("NO_SOURCE_EXPORT", self.payload["status"])
        self.assertIn("NO_BRIDGE", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_failure_cover_counts(self) -> None:
        self.assertTrue(self.theorem["p2534_boolean_prime_implicant_inherited"])
        self.assertEqual(self.theorem["valid_ternary_assignment_count"], 81)
        self.assertEqual(self.theorem["strict_accept_assignment_count"], 1)
        self.assertEqual(self.theorem["failure_assignment_count"], 80)
        self.assertEqual(self.theorem["strict_theorem_failure_cover_term_count"], 4)
        self.assertTrue(self.theorem["each_strict_theorem_cover_term_covers_54_rows"])
        self.assertTrue(self.theorem["each_strict_theorem_cover_term_splits_27_axiom_27_absent"])
        for row in self.cert["strict_theorem_failure_cover_rows"]:
            self.assertEqual(row["covered_failure_rows"], 54)
            self.assertEqual(row["covered_axiom_rows"], 27)
            self.assertEqual(row["covered_absent_rows"], 27)

    def test_intersections_and_union(self) -> None:
        self.assertTrue(self.theorem["cover_intersections_match_closed_form"])
        by_size = {}
        for row in self.cert["strict_theorem_failure_cover_intersections"]:
            by_size.setdefault(row["subset_size"], set()).add(row["intersection_failure_rows"])
        self.assertEqual(by_size, {1: {54}, 2: {36}, 3: {24}, 4: {16}})
        self.assertEqual(self.theorem["strict_theorem_failure_cover_union_count_by_inclusion_exclusion"], 80)
        self.assertTrue(self.theorem["strict_theorem_failure_cover_exhausts_all_failures"])

    def test_prime_failure_implicants_and_nearest_rows(self) -> None:
        self.assertEqual(self.theorem["exact_prime_failure_implicant_count"], 8)
        self.assertTrue(self.theorem["strict_missing_prime_cubes_present"])
        self.assertTrue(self.theorem["absent_prime_cubes_present_due_valid_embedding"])
        self.assertEqual(self.theorem["nearest_failure_row_count"], 8)
        prime_cubes = {row["cube"] for row in self.cert["exact_prime_failure_implicants"]}
        self.assertIn("t_M_strict=0", prime_cubes)
        self.assertIn("p_M_present=0", prime_cubes)
        self.assertEqual(len(self.cert["nearest_failure_rows"]), 8)

    def test_negative_controls(self) -> None:
        self.assertTrue(self.theorem["dual_failure_cover_certificate_exported"])
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
        self.assertIn("P2535/S1485", MD.read_text(encoding="utf-8"))
        self.assertIn("P2535/S1485", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2535/S1485", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
