from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2692_s1642_target_independent_positive_beta_zbeta_source_audit.py"
OUT = ROOT / "generated" / "p2692_s1642_target_independent_positive_beta_zbeta_source_audit.json"
MD = ROOT / "generated" / "p2692_s1642_target_independent_positive_beta_zbeta_source_audit.md"


class P2692TargetIndependentPositiveBetaZBetaSourceAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_grep_covers_beta_zbeta_frontier(self) -> None:
        audit = self.payload["content_grep"]
        self.assertIn("positive beta/Z_beta", audit["mode"])
        for key in (
            "p2691_selected_p2692",
            "beta_zbeta_obstructions",
            "normalization_orbit_contract",
            "typed_metric_uv_obligations",
            "forbidden_imports",
        ):
            self.assertIn(key, audit["patterns"])
            self.assertGreater(audit["patterns"][key]["count"], 0)

    def test_state_reads_keep_beta_source_unexported(self) -> None:
        state = self.payload["state_reads"]
        self.assertTrue(state["p2691_selected_p2692"])
        self.assertFalse(state["p2680_positive_beta_zbeta_source_exported"])
        self.assertFalse(state["p2680_canonical_length_uv_source_exported"])
        self.assertFalse(state["p2649_beta_source_exported"])
        self.assertEqual(state["p2649_passing_routes"], [])
        self.assertTrue(state["p2651_beta_one_working_gauge_allowed"])
        self.assertFalse(state["p2651_beta_source_exported"])
        self.assertTrue(state["p2653_scale_orbit_equivalence_verified"])
        self.assertFalse(state["p2653_beta_source_exported"])
        self.assertFalse(state["p2653_canonical_unit_exported"])

    def test_beta_orbit_witness_is_positive_but_not_source(self) -> None:
        orbit = self.payload["beta_orbit_witness"]
        self.assertTrue(orbit["all_positive_betas_have_beta_one_representative"])
        self.assertLess(orbit["max_abs_invariant_residual"], 1e-10)
        self.assertFalse(orbit["source_generated_by_orbit"])
        self.assertEqual(len(orbit["rows"]), 6)

    def test_tail_ratio_inversion_is_target_dependent(self) -> None:
        witness = self.payload["tail_ratio_inversion_witness"]
        self.assertTrue(witness["positive_beta_recoverable_for_multiple_declared_targets"])
        self.assertFalse(witness["target_independent_source_generated"])
        self.assertEqual(len(witness["rows"]), 4)

    def test_candidate_matrix_has_no_passing_source(self) -> None:
        matrix = self.payload["source_candidate_matrix"]
        self.assertEqual(len(matrix), 5)
        self.assertFalse(any(row["positive_beta_available"] and row["target_independent"] and row["source_exported_now"] for row in matrix))
        self.assertTrue(any(row["candidate"] == "normalization_orbit_beta_equals_one_representative" and row["positive_beta_available"] and not row["target_independent"] for row in matrix))
        self.assertTrue(any(row["candidate"] == "dimensionless_conservation_or_operator_identity_unique_beta_one" and row["target_independent"] and not row["source_exported_now"] for row in matrix))

    def test_decision_freezes_beta_atom_and_updates_docs(self) -> None:
        decision = self.payload["decision"]
        self.assertFalse(decision["target_independent_positive_beta_or_zbeta_source_exported_now"])
        self.assertTrue(decision["bounded_no_go_for_current_beta_zbeta_atom"])
        self.assertIn("P2693", decision["next_honest_step"])
        self.assertFalse(decision["role_transfer_started_now"])
        self.assertFalse(decision["ltotal_promoted_now"])
        self.assertFalse(decision["toe_closed_now"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2692/S1642", MD.read_text(encoding="utf-8"))
        self.assertIn("P2692/S1642", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2692/S1642", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2692/S1642", (ROOT.parent / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
