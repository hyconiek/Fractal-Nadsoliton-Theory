from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2551_s1501_strict_damping_post_prime_log_slope_anchor_obstruction_certificate import (
    MD,
    OUT,
    PRIMES,
    append_doc_sections,
    build_payload,
    write_markdown,
)


class P2551StrictDampingPostPrimeLogSlopeAnchorObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["strict_damping_post_prime_log_slope_anchor_obstruction_certificate"]["theorem_export"]
        cls.audit = cls.theorem["post_prime_log_slope_anchor_obstruction_audit"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2551")
        self.assertEqual(self.payload["stage_id"], "S1501")
        self.assertIn("POST_PRIME_LOG_SLOPE_ANCHOR_OBSTRUCTION", self.payload["status"])
        self.assertIn("NO_SLOPE_SOURCE_EXPORT", self.payload["status"])
        self.assertIn("NO_BRIDGE", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_precursors_inherited(self) -> None:
        self.assertTrue(self.theorem["p2528_prime_slope_anchor_equivalence_inherited"])
        self.assertTrue(self.theorem["p2543_slope_value_current_premise_obstruction_inherited"])
        self.assertTrue(self.theorem["p2547_post_identity_residual_trikey_inherited"])
        self.assertTrue(self.theorem["p2550_prime_log_adjacent_basis_inherited"])

    def test_post_prime_log_slope_countermodels(self) -> None:
        rows = self.theorem["post_prime_log_slope_line_rows"]
        self.assertEqual(len(rows), 7)
        self.assertTrue(all(row["satisfies_p2550_adjacent_ratio_basis"] for row in rows))
        self.assertEqual(self.audit["strict_accepting_slopes"], ["4/5"])
        self.assertEqual(self.audit["countermodel_count"], 6)
        self.assertIn("1/2", self.audit["countermodel_slopes"])
        countermodel = self.audit["explicit_delta_one_half_countermodel"]
        self.assertEqual(countermodel["slope_delta"], "1/2")
        self.assertTrue(countermodel["satisfies_p2550_adjacent_ratio_basis"])
        self.assertFalse(countermodel["strict_delta_4_over_5_anchor_accepts"])
        self.assertTrue(countermodel["post_prime_log_slope_countermodel"])
        self.assertFalse(self.audit["prime_log_adjacent_basis_entails_slope_value"])
        self.assertTrue(self.theorem["post_prime_log_slope_value_nonentailment_exported"])

    def test_single_prime_anchor_equivalence(self) -> None:
        anchors = self.audit["single_prime_anchor_rows"]
        self.assertEqual(len(anchors), len(PRIMES))
        self.assertEqual([row["anchor_prime"] for row in anchors], PRIMES)
        self.assertTrue(all(row["determinant_log_prime_nonzero"] for row in anchors))
        self.assertTrue(all(row["single_prime_anchor_suffices_given_prime_log_line"] for row in anchors))
        self.assertTrue(self.audit["all_single_prime_anchors_suffice_given_prime_log_line"])
        self.assertTrue(self.theorem["single_prime_anchor_equivalence_reconfirmed_given_prime_log_line"])

    def test_negative_controls_and_recommendation(self) -> None:
        self.assertFalse(self.theorem["slope_value_or_prime_anchor_source_exported"])
        self.assertFalse(self.theorem["prime_log_proportionality_source_exported"])
        self.assertFalse(self.theorem["m2_operator_signature_source_exported"])
        self.assertFalse(self.theorem["strict_quadruple_trace_source_exported"])
        self.assertFalse(self.theorem["beta_eta_numeric_source_exported"])
        self.assertFalse(self.theorem["strict_damping_beta_eta_source_exported"])
        self.assertFalse(self.theorem["source_obligation_discharge_exported"])
        self.assertFalse(self.theorem["full_bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["selector_closure_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_certificate"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertTrue(self.theorem["conditional_full_strict_damping_if_all_other_sources_and_slope_anchor_supplied"])
        self.assertEqual(self.theorem["residual_after_hypothetical_identity_m2_prime_log_and_slope"], [])
        self.assertIn("prime-value/slope-anchor theorem", self.theorem["recommended_next_honest_step"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

    def test_rg_audit_and_documentation_markers(self) -> None:
        audit = self.payload["rg_non_duplication_audit"]
        self.assertEqual(audit["tool"], "rg")
        self.assertIn("intended_research_nonduplication", audit["patterns"])
        self.assertGreaterEqual(audit["patterns"]["new_packet"]["count"], 0)
        self.assertIn("P2551/S1501", MD.read_text(encoding="utf-8"))
        self.assertIn("P2551/S1501", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2551/S1501", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
