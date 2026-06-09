from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2602_s1552_nadsoliton_rg_fixed_rate_prime_log_source_theorem import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2602NadsolitonRGFixedRatePrimeLogSourceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["nadsoliton_rg_fixed_rate_prime_log_source_theorem"]["theorem_export"]
        cls.derivation = cls.theorem["rg_fixed_rate_prime_log_derivation"]
        cls.residual = cls.theorem["post_prime_log_post_unital_post_m2_residual_matrix"]

    def test_identity_and_inheritance(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2602")
        self.assertEqual(self.payload["stage_id"], "S1552")
        self.assertIn("PRIME_LOG_SOURCE", self.payload["status"])
        self.assertTrue(self.payload["gatekeeper_checks"]["p2542_obstruction_inherited"])
        self.assertTrue(self.payload["gatekeeper_checks"]["p2601_unital_multiplicative_source_inherited"])

    def test_fixed_rate_sources_prime_log_but_not_slope(self) -> None:
        self.assertTrue(self.theorem["rg_fixed_rate_prime_log_source_exported"])
        self.assertTrue(self.theorem["prime_log_proportionality_source_exported"])
        self.assertTrue(self.derivation["all_gamma_samples_prime_log_proportional"])
        self.assertTrue(self.derivation["non_strict_gamma_samples_also_prime_log_proportional"])
        self.assertTrue(self.derivation["strict_slope_not_selected_by_fixed_rate_alone"])
        self.assertEqual(self.derivation["prime_generator_formula"], "v_p = y_p = gamma_star log(p)")
        self.assertEqual(len(self.derivation["gamma_sample_rows"]), 3)
        for row in self.derivation["gamma_sample_rows"]:
            self.assertLess(row["ratio_spread"], 1e-12)
            self.assertLess(row["max_multiplicative_character_defect"], 1e-12)

    def test_post_prime_log_residual_matrix(self) -> None:
        self.assertEqual(self.residual["remaining_keys_after_m2_unital_prime_log"], ["slope_value_or_prime_anchor_source"])
        self.assertEqual(self.residual["residual_truth_table_row_count"], 2)
        self.assertEqual(self.residual["residual_accepting_row_count"], 1)
        self.assertEqual(self.residual["remaining_missing_source_key_count_after_p2602"], 1)
        self.assertFalse(self.residual["strict_damping_beta_eta_source_accepts_after_p2602_current_assignment"])
        self.assertFalse(self.theorem["strict_damping_beta_eta_source_exported"])

    def test_scope_guards_and_docs(self) -> None:
        self.assertFalse(self.theorem["slope_value_or_prime_anchor_source_exported"])
        self.assertFalse(self.theorem["beta_eta_numeric_source_exported"])
        self.assertFalse(self.theorem["bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_theorem"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2602/S1552", MD.read_text(encoding="utf-8"))
        self.assertIn("P2602/S1552", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2602/S1552", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
