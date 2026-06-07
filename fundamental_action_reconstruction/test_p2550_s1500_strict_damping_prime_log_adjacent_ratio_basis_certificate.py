from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2550_s1500_strict_damping_prime_log_adjacent_ratio_basis_certificate import (
    MD,
    OUT,
    PRIMES,
    append_doc_sections,
    build_payload,
    write_markdown,
)


class P2550StrictDampingPrimeLogAdjacentRatioBasisTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["strict_damping_prime_log_adjacent_ratio_basis_certificate"]["theorem_export"]
        cls.basis = cls.theorem["prime_log_adjacent_ratio_basis_audit"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2550")
        self.assertEqual(self.payload["stage_id"], "S1500")
        self.assertIn("PRIME_LOG_ADJACENT_RATIO_BASIS", self.payload["status"])
        self.assertIn("SOURCE_OBLIGATION_ONLY", self.payload["status"])
        self.assertIn("NO_PRIME_LOG_SOURCE_EXPORT", self.payload["status"])
        self.assertIn("NO_BRIDGE", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_precursors_inherited(self) -> None:
        self.assertTrue(self.theorem["p2527_prime_log_slope_line_inherited"])
        self.assertTrue(self.theorem["p2542_prime_log_current_premise_obstruction_inherited"])
        self.assertTrue(self.theorem["p2547_post_identity_residual_trikey_inherited"])
        self.assertTrue(self.theorem["p2549_trace_source_obstruction_inherited"])

    def test_adjacent_constraint_basis(self) -> None:
        self.assertEqual(self.basis["prime_order"], PRIMES)
        self.assertEqual(self.basis["adjacent_edges_by_prime"], [[2, 3], [3, 5], [5, 7], [7, 11]])
        self.assertEqual(self.basis["constraint_matrix_rows"], [
            [1, -1, 0, 0, 0],
            [0, 1, -1, 0, 0],
            [0, 0, 1, -1, 0],
            [0, 0, 0, 1, -1],
        ])
        self.assertEqual(self.basis["constraint_matrix_rank"], 4)
        self.assertEqual(self.basis["constraint_matrix_nullity"], 1)
        self.assertEqual(self.basis["nullspace_basis"], [[1, 1, 1, 1, 1]])
        self.assertTrue(self.basis["nullspace_is_constant_ratio_line"])
        self.assertTrue(self.basis["full_adjacent_basis_equivalent_to_prime_log_proportionality"])
        self.assertEqual(self.theorem["minimal_constraint_count_for_five_prime_ratios"], 4)
        self.assertTrue(self.theorem["full_adjacent_ratio_basis_suffices_for_prime_log_proportionality"])

    def test_single_omission_witnesses(self) -> None:
        audit = self.basis["proper_subset_audit"]
        self.assertEqual(audit["single_omission_witness_count"], 4)
        self.assertTrue(audit["all_single_omissions_have_nonproportional_witness"])
        self.assertTrue(audit["all_single_omissions_satisfy_remaining_edges"])
        self.assertTrue(audit["basis_irredundant_under_single_omission"])
        for witness in audit["single_omission_witnesses"]:
            self.assertTrue(witness["satisfies_all_non_omitted_adjacent_equalities"])
            self.assertTrue(witness["violates_omitted_adjacent_equality"])
            self.assertEqual(witness["ratio_spread"], 1)
            self.assertFalse(witness["prime_log_proportionality_accepts"])
        self.assertTrue(self.theorem["single_omission_countermodels_exported"])

    def test_negative_controls_and_recommendation(self) -> None:
        self.assertTrue(self.theorem["adjacent_ratio_basis_is_source_obligation_not_source"])
        self.assertFalse(self.theorem["prime_log_proportionality_source_exported"])
        self.assertFalse(self.theorem["slope_value_or_prime_anchor_source_exported"])
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
        self.assertEqual(self.theorem["residual_after_hypothetical_identity_m2_and_prime_log"], ["slope_value_or_prime_anchor_source"])
        self.assertIn("slope_value_or_prime_anchor_source", self.theorem["recommended_next_honest_step"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

    def test_rg_audit_and_documentation_markers(self) -> None:
        audit = self.payload["rg_non_duplication_audit"]
        self.assertEqual(audit["tool"], "rg")
        self.assertIn("intended_research_nonduplication", audit["patterns"])
        self.assertGreaterEqual(audit["patterns"]["new_packet"]["count"], 0)
        self.assertIn("P2550/S1500", MD.read_text(encoding="utf-8"))
        self.assertIn("P2550/S1500", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2550/S1500", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
