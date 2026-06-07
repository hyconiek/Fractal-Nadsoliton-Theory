from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2548_s1498_strict_damping_m2_trace_arity_conditional_selection_certificate import (
    MD,
    OUT,
    TARGET_TRACE_ARITY,
    append_doc_sections,
    build_payload,
    write_markdown,
)


class P2548StrictDampingM2TraceArityConditionalSelectionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["strict_damping_m2_trace_arity_conditional_selection_certificate"]["theorem_export"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2548")
        self.assertEqual(self.payload["stage_id"], "S1498")
        self.assertIn("M2_TRACE_ARITY_CONDITIONAL_SELECTION", self.payload["status"])
        self.assertIn("NO_ACTUAL_TRACE_SOURCE_EXPORT", self.payload["status"])
        self.assertIn("NO_BRIDGE", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_precursors_inherited(self) -> None:
        self.assertTrue(self.theorem["p2540_m2_obstruction_inherited"])
        self.assertTrue(self.theorem["p2544_four_key_blocker_inherited"])
        self.assertTrue(self.theorem["p2547_post_identity_residual_trikey_inherited"])

    def test_trace_arity_selection(self) -> None:
        selection = self.theorem["trace_arity_selection_audit"]
        self.assertEqual(selection["target_trace_arity"], TARGET_TRACE_ARITY)
        self.assertEqual(selection["selection_equation"], "2*m = 4")
        self.assertEqual(selection["integer_solution_m"], 2)
        self.assertTrue(selection["equation_has_integer_solution"])
        self.assertEqual(selection["matching_orders_in_audited_range"], [2])
        self.assertTrue(selection["unique_matching_order_in_audited_range"])
        self.assertEqual(selection["nonmatching_order_witness_count"], 9)
        rows = self.theorem["trace_arity_rows"]
        self.assertEqual(len(rows), 10)
        self.assertTrue(all(row["self_adjoint_boundary_trace_arity"] == 2 * row["derivative_order_m"] for row in rows))
        self.assertEqual([row["derivative_order_m"] for row in rows if row["matches_strict_quadruple_trace_arity"]], [2])
        self.assertTrue(self.theorem["conditional_quadruple_trace_selects_m2"])
        self.assertTrue(self.theorem["conditional_m2_closure_if_quadruple_trace_source_supplied"])

    def test_conditional_not_actual_source(self) -> None:
        self.assertTrue(self.theorem["strict_quadruple_trace_source_required"])
        self.assertFalse(self.theorem["strict_quadruple_trace_source_exported"])
        self.assertFalse(self.theorem["m2_operator_signature_source_exported"])
        self.assertIn("does not prove", self.theorem["conditional_not_actual_source_reason"])

    def test_residual_after_identity_and_conditional_m2(self) -> None:
        self.assertEqual(self.theorem["residual_keys_after_identity_and_conditional_m2"], [
            "prime_log_proportionality_source",
            "slope_value_or_prime_anchor_source",
        ])
        self.assertEqual(self.theorem["residual_truth_table_after_identity_and_conditional_m2_row_count"], 4)
        self.assertEqual(self.theorem["residual_accepting_row_count_after_identity_and_conditional_m2"], 1)
        accepting = self.theorem["residual_accepting_row_after_identity_and_conditional_m2"]
        self.assertTrue(all(accepting["residual_assignment"].values()))
        self.assertTrue(accepting["strict_damping_beta_eta_source_accepts"])
        self.assertTrue(self.theorem["conditional_identity_and_m2_reduce_missing_count_from_3_to_2"])
        self.assertTrue(self.theorem["identity_plus_conditional_m2_still_cannot_export_beta_eta_numeric_source"])
        self.assertTrue(self.theorem["identity_plus_conditional_m2_still_cannot_export_strict_damping_beta_eta_source"])

    def test_negative_controls_and_recommendation(self) -> None:
        self.assertFalse(self.theorem["multiplicative_character_law_source_exported"])
        self.assertFalse(self.theorem["prime_log_proportionality_source_exported"])
        self.assertFalse(self.theorem["slope_value_or_prime_anchor_source_exported"])
        self.assertFalse(self.theorem["beta_eta_numeric_source_exported"])
        self.assertFalse(self.theorem["strict_damping_beta_eta_source_exported"])
        self.assertFalse(self.theorem["source_obligation_discharge_exported"])
        self.assertFalse(self.theorem["full_bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["selector_closure_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_certificate"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertIn("quadruple-boundary-trace theorem", self.theorem["recommended_next_honest_step"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

    def test_rg_audit_and_documentation_markers(self) -> None:
        audit = self.payload["rg_non_duplication_audit"]
        self.assertEqual(audit["tool"], "rg")
        self.assertIn("intended_research_nonduplication", audit["patterns"])
        self.assertGreaterEqual(audit["patterns"]["new_packet"]["count"], 0)
        self.assertIn("P2548/S1498", MD.read_text(encoding="utf-8"))
        self.assertIn("P2548/S1498", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2548/S1498", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
