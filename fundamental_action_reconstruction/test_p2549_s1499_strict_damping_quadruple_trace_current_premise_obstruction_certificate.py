from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2549_s1499_strict_damping_quadruple_trace_current_premise_obstruction_certificate import (
    MD,
    OUT,
    TARGET_TRACE_ARITY,
    append_doc_sections,
    build_payload,
    write_markdown,
)


class P2549StrictDampingQuadrupleTraceCurrentPremiseObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["strict_damping_quadruple_trace_current_premise_obstruction_certificate"]["theorem_export"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2549")
        self.assertEqual(self.payload["stage_id"], "S1499")
        self.assertIn("QUADRUPLE_TRACE_CURRENT_PREMISE_OBSTRUCTION", self.payload["status"])
        self.assertIn("NO_TRACE_SOURCE_EXPORT", self.payload["status"])
        self.assertIn("NO_M2_SOURCE_EXPORT", self.payload["status"])
        self.assertIn("NO_BRIDGE", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_precursors_inherited(self) -> None:
        self.assertTrue(self.theorem["p2540_current_premise_order_nonidentifiability_inherited"])
        self.assertTrue(self.theorem["p2547_post_identity_residual_trikey_inherited"])
        self.assertTrue(self.theorem["p2548_conditional_trace_arity_selection_inherited"])
        self.assertTrue(self.theorem["p2548_actual_trace_source_absent_inherited"])

    def test_quadruple_trace_nonentailment_rows(self) -> None:
        rows = self.theorem["trace_nonentailment_rows"]
        self.assertEqual(len(rows), 10)
        self.assertTrue(all(row["inherits_current_source_free_premise_pass"] for row in rows))
        self.assertEqual([row["self_adjoint_trace_arity"] for row in rows], list(range(2, 21, 2)))
        self.assertEqual([row["derivative_order_m"] for row in rows if row["satisfies_exact_quadruple_trace"]], [2])
        self.assertEqual([row["derivative_order_m"] for row in rows if row["countermodel_to_quadruple_trace_source"]], [1, 3, 4, 5, 6, 7, 8, 9, 10])

    def test_obstruction_audit(self) -> None:
        audit = self.theorem["trace_obstruction_audit"]
        self.assertEqual(audit["target_trace_arity"], TARGET_TRACE_ARITY)
        self.assertEqual(audit["current_premise_passing_order_count"], 10)
        self.assertEqual(audit["quadruple_trace_accepting_orders"], [2])
        self.assertEqual(audit["quadruple_trace_accepting_order_count"], 1)
        self.assertEqual(audit["quadruple_trace_countermodel_count"], 9)
        self.assertEqual(audit["m3_countermodel"]["derivative_order_m"], 3)
        self.assertEqual(audit["m3_countermodel"]["self_adjoint_trace_arity"], 6)
        self.assertTrue(audit["m3_countermodel"]["countermodel_to_quadruple_trace_source"])
        self.assertFalse(audit["current_premises_entail_exact_quadruple_trace"])
        self.assertFalse(audit["current_premises_entail_m2_via_quadruple_trace"])
        self.assertTrue(self.theorem["current_premise_quadruple_trace_nonentailment_exported"])
        self.assertTrue(self.theorem["current_premise_m2_via_quadruple_trace_route_refuted"])

    def test_negative_controls_and_recommendation(self) -> None:
        self.assertTrue(self.theorem["strict_quadruple_trace_source_required_for_p2548"])
        self.assertFalse(self.theorem["strict_quadruple_trace_source_exported"])
        self.assertFalse(self.theorem["m2_operator_signature_source_exported"])
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
        self.assertEqual(self.theorem["residual_after_failed_trace_source"], [
            "strict_quadruple_trace_source",
            "prime_log_proportionality_source",
            "slope_value_or_prime_anchor_source",
        ])
        self.assertIn("prime_log_proportionality_source", self.theorem["recommended_next_honest_step"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

    def test_rg_audit_and_documentation_markers(self) -> None:
        audit = self.payload["rg_non_duplication_audit"]
        self.assertEqual(audit["tool"], "rg")
        self.assertIn("intended_research_nonduplication", audit["patterns"])
        self.assertGreaterEqual(audit["patterns"]["new_packet"]["count"], 0)
        self.assertIn("P2549/S1499", MD.read_text(encoding="utf-8"))
        self.assertIn("P2549/S1499", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2549/S1499", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
