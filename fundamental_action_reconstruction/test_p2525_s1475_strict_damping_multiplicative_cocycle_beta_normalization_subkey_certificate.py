from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2525_s1475_strict_damping_multiplicative_cocycle_beta_normalization_subkey_certificate import (
    MD,
    OUT,
    build_payload,
    write_markdown,
)


class P2525StrictDampingMultiplicativeCocycleBetaNormalizationSubkeyTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["strict_damping_multiplicative_cocycle_beta_normalization_subkey_certificate"]["theorem_export"]
        cls.cert = cls.theorem["strict_damping_multiplicative_cocycle_beta_normalization_subkey_certificate"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2525")
        self.assertEqual(self.payload["stage_id"], "S1475")
        self.assertIn("MULTIPLICATIVE_COCYCLE_BETA_NORMALIZATION_SUBKEY", self.payload["status"])
        self.assertIn("CONDITIONAL_SUBKEY_ONLY", self.payload["status"])
        self.assertIn("NO_SLOPE_SOURCE", self.payload["status"])
        self.assertIn("NO_OPERATOR_SOURCE", self.payload["status"])
        self.assertIn("NO_BRIDGE", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_multiplicative_cocycle_pins_beta_normalization_only(self) -> None:
        self.assertTrue(self.theorem["p2524_continuum_nonidentifiability_inherited"])
        self.assertEqual(self.theorem["candidate_grid_row_count"], 35)
        self.assertGreater(self.theorem["multiplicative_pair_count_on_domain_1_to_11"], 0)
        self.assertEqual(self.theorem["multiplicative_character_accepting_row_count"], 7)
        self.assertTrue(self.theorem["multiplicative_law_pins_log_beta_zero_on_grid"])
        self.assertTrue(self.theorem["multiplicative_law_leaves_slope_continuum_on_grid"])
        self.assertTrue(self.theorem["conditional_beta_normalization_subkey_exported"])
        accepted = [row for row in self.cert["multiplicative_cocycle_rows"] if row["multiplicative_character_accepts"]]
        self.assertEqual({row["intercept_log_beta"] for row in accepted}, {"0"})
        self.assertEqual(len({row["slope_delta"] for row in accepted}), 7)
        self.assertTrue(all(row["defect_equals_minus_intercept_on_all_pairs"] for row in self.cert["multiplicative_cocycle_rows"]))

    def test_strict_target_still_needs_slope_filter(self) -> None:
        audit = self.cert["multiplicative_filter_audit"]
        self.assertTrue(audit["multiplicative_law_pins_log_beta_zero_on_grid"])
        self.assertTrue(audit["multiplicative_law_leaves_slope_continuum_on_grid"])
        self.assertEqual(audit["strict_slope_filter_after_multiplicative_law_count"], 1)
        self.assertTrue(audit["strict_slope_filter_after_multiplicative_law_unique_target"])
        self.assertIn("slope_value_source", audit["source_obligation_normal_form"])

    def test_negative_controls(self) -> None:
        self.assertTrue(self.theorem["conditional_beta_normalization_subkey_exported"])
        self.assertFalse(self.theorem["multiplicative_character_law_source_exported"])
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
        self.assertIn("P2525/S1475", MD.read_text(encoding="utf-8"))
        self.assertIn("P2525/S1475", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2525/S1475", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
