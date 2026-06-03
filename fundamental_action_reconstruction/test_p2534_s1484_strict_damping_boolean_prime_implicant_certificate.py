from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2534_s1484_strict_damping_boolean_prime_implicant_certificate import (
    MD,
    OUT,
    build_payload,
    write_markdown,
)


class P2534StrictDampingBooleanPrimeImplicantTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["strict_damping_boolean_prime_implicant_certificate"]["theorem_export"]
        cls.cert = cls.theorem["strict_damping_boolean_prime_implicant_certificate"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2534")
        self.assertEqual(self.payload["stage_id"], "S1484")
        self.assertIn("BOOLEAN_PRIME_IMPLICANT", self.payload["status"])
        self.assertIn("UNIQUE_ALL_STRICT_CUBE", self.payload["status"])
        self.assertIn("NO_SOURCE_EXPORT", self.payload["status"])
        self.assertIn("NO_BRIDGE", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_boolean_embedding_counts(self) -> None:
        self.assertTrue(self.theorem["p2533_polynomial_certificate_inherited"])
        self.assertEqual(self.theorem["valid_ternary_assignment_count"], 81)
        self.assertEqual(self.theorem["all_boolean_cube_count"], 6561)
        self.assertEqual(self.theorem["strict_true_assignment_count"], 1)
        self.assertEqual(self.theorem["strict_implicant_count"], 16)
        self.assertEqual(self.cert["strict_implicant_degree_histogram"], {
            "0": 0,
            "1": 0,
            "2": 0,
            "3": 0,
            "4": 1,
            "5": 4,
            "6": 6,
            "7": 4,
            "8": 1,
        })

    def test_unique_prime_implicant(self) -> None:
        self.assertEqual(self.theorem["minimal_strict_implicant_degree"], 4)
        self.assertTrue(self.theorem["minimal_strict_implicant_degree_is_four"])
        self.assertEqual(len(self.theorem["minimal_strict_implicants"]), 1)
        self.assertTrue(self.theorem["unique_prime_strict_implicant"])
        self.assertTrue(self.theorem["prime_implicant_is_all_four_strict_theorems"])
        prime = self.theorem["prime_strict_implicants"][0]
        self.assertEqual(prime["degree"], 4)
        self.assertEqual(set(prime["literals"]), {"t_M_strict", "t_P_strict", "t_A_strict", "t_O_strict"})
        self.assertEqual(set(prime["literals"].values()), {1})
        self.assertEqual(self.theorem["p_only_strict_implicant_count"], 0)
        self.assertTrue(self.theorem["no_present_only_cube_implies_strict_acceptance"])

    def test_essential_literal_witnesses(self) -> None:
        self.assertEqual(self.theorem["omitted_strict_literal_witness_count"], 4)
        self.assertTrue(self.theorem["every_strict_literal_has_axiom_and_absent_false_witnesses"])
        self.assertEqual(self.theorem["nearest_axiom_false_rows_count"], 4)
        omitted = {row["omitted_literal"] for row in self.cert["omitted_strict_literal_witnesses"]}
        self.assertEqual(omitted, {"t_M_strict", "t_P_strict", "t_A_strict", "t_O_strict"})
        for row in self.cert["omitted_strict_literal_witnesses"]:
            self.assertEqual(row["false_rows_after_relaxation_count"], 2)
            self.assertIn("axiom", set(row["axiom_augmented_false_witness"].values()))
            self.assertIn("absent", set(row["blocked_missing_key_false_witness"].values()))

    def test_negative_controls(self) -> None:
        self.assertTrue(self.theorem["boolean_prime_implicant_certificate_exported"])
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
        self.assertIn("P2534/S1484", MD.read_text(encoding="utf-8"))
        self.assertIn("P2534/S1484", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2534/S1484", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
