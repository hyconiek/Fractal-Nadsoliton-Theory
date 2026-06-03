from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2526_s1476_strict_damping_finite_monoid_prime_character_nullity_certificate import (
    MD,
    OUT,
    build_payload,
    write_markdown,
)


class P2526StrictDampingFiniteMonoidPrimeCharacterNullityTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["strict_damping_finite_monoid_prime_character_nullity_certificate"]["theorem_export"]
        cls.cert = cls.theorem["strict_damping_finite_monoid_prime_character_nullity_certificate"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2526")
        self.assertEqual(self.payload["stage_id"], "S1476")
        self.assertIn("FINITE_MONOID_PRIME_CHARACTER_NULLITY", self.payload["status"])
        self.assertIn("CONDITIONAL_CHARACTER_FAMILY_ONLY", self.payload["status"])
        self.assertIn("NO_PRIME_PROPORTIONALITY_SOURCE", self.payload["status"])
        self.assertIn("NO_SLOPE_SOURCE", self.payload["status"])
        self.assertIn("NO_OPERATOR_SOURCE", self.payload["status"])
        self.assertIn("NO_BRIDGE", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_exact_rank_nullity_and_prime_free_variables(self) -> None:
        self.assertTrue(self.theorem["p2525_beta_normalization_subkey_inherited"])
        self.assertEqual(self.theorem["multiplicative_pair_count_on_domain_1_to_11"], 29)
        self.assertEqual(self.theorem["constraint_row_count"], 29)
        self.assertEqual(self.theorem["constraint_column_count"], 11)
        self.assertEqual(self.theorem["exact_constraint_rank"], 6)
        self.assertEqual(self.theorem["exact_constraint_nullity"], 5)
        self.assertEqual(self.theorem["canonical_prime_parameter_variables_y_indices_1_based"], [2, 3, 5, 7, 11])
        self.assertTrue(self.theorem["factorization_parameterization_dimension_matches_nullity"])
        self.assertTrue(self.theorem["multiplicative_law_alone_leaves_prime_generator_freedom"])

    def test_multiplicativity_is_not_affine_slope_source(self) -> None:
        rows = self.cert["sample_prime_character_rows"]
        self.assertTrue(self.theorem["strict_log_character_is_one_member_of_prime_character_family"])
        self.assertTrue(self.theorem["non_affine_prime_characters_pass_multiplicativity"])
        self.assertTrue(all(row["multiplicative_constraints_accept"] for row in rows))
        non_affine = [row for row in rows if not row["affine_secants_constant"]]
        self.assertGreaterEqual(len(non_affine), 1)
        self.assertTrue(all(row["secant_spread"] > 0 for row in non_affine))
        self.assertTrue(self.theorem["prime_log_proportionality_needed_to_recover_affine_slope_line"])
        self.assertTrue(self.theorem["slope_value_source_still_needed_for_delta_4_over_5"])

    def test_negative_controls(self) -> None:
        self.assertTrue(self.theorem["conditional_finite_monoid_character_nullity_exported"])
        self.assertFalse(self.theorem["multiplicative_character_law_source_exported"])
        self.assertFalse(self.theorem["prime_log_proportionality_source_exported"])
        self.assertFalse(self.theorem["slope_value_source_exported"])
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
        self.assertIn("P2526/S1476", MD.read_text(encoding="utf-8"))
        self.assertIn("P2526/S1476", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2526/S1476", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
