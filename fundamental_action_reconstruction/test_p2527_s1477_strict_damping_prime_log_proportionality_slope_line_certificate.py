from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2527_s1477_strict_damping_prime_log_proportionality_slope_line_certificate import (
    MD,
    OUT,
    build_payload,
    write_markdown,
)


class P2527StrictDampingPrimeLogProportionalitySlopeLineTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["strict_damping_prime_log_proportionality_slope_line_certificate"]["theorem_export"]
        cls.cert = cls.theorem["strict_damping_prime_log_proportionality_slope_line_certificate"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2527")
        self.assertEqual(self.payload["stage_id"], "S1477")
        self.assertIn("PRIME_LOG_PROPORTIONALITY_SLOPE_LINE", self.payload["status"])
        self.assertIn("CONDITIONAL_SUBKEY_ONLY", self.payload["status"])
        self.assertIn("NO_SLOPE_VALUE_SOURCE", self.payload["status"])
        self.assertIn("NO_OPERATOR_SOURCE", self.payload["status"])
        self.assertIn("NO_BRIDGE", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_prime_log_proportionality_rank_nullity(self) -> None:
        self.assertTrue(self.theorem["p2526_prime_character_nullity_inherited"])
        self.assertEqual(self.theorem["exact_normalized_ratio_constraint_rank"], 4)
        self.assertEqual(self.theorem["exact_normalized_ratio_constraint_nullity"], 1)
        self.assertTrue(self.theorem["one_dimensional_slope_line_after_prime_log_proportionality"])
        self.assertEqual(self.theorem["slope_candidate_count"], 7)
        self.assertTrue(self.theorem["all_slope_candidates_pass_prime_log_proportionality"])
        self.assertTrue(self.theorem["all_slope_candidates_give_affine_node_line"])

    def test_strict_target_is_not_unique_without_slope_source(self) -> None:
        rows = self.cert["slope_line_rows"]
        strict_rows = [row for row in rows if row["is_strict_numeric_target"]]
        self.assertEqual(len(strict_rows), 1)
        self.assertTrue(self.theorem["strict_numeric_target_is_one_slope_line_member"])
        self.assertGreater(len(rows), 1)
        self.assertTrue(all(row["prime_log_proportionality_accepts"] for row in rows))
        self.assertTrue(all(row["affine_node_line_accepts"] for row in rows))

    def test_nonproportional_prime_characters_rejected(self) -> None:
        rejected = self.cert["representative_nonproportional_rejection_rows"]
        self.assertTrue(self.theorem["representative_nonproportional_prime_characters_rejected"])
        self.assertGreaterEqual(len(rejected), 3)
        self.assertTrue(all(row["prime_log_proportionality_rejects"] for row in rejected))
        self.assertTrue(all(row["affine_node_line_rejects"] for row in rejected))

    def test_negative_controls(self) -> None:
        self.assertTrue(self.theorem["prime_log_proportionality_subkey_exported"])
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
        self.assertIn("P2527/S1477", MD.read_text(encoding="utf-8"))
        self.assertIn("P2527/S1477", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2527/S1477", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
