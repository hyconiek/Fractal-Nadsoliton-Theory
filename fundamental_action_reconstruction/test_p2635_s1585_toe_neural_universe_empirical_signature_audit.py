from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2635_s1585_toe_neural_universe_empirical_signature_audit.py"
OUT = ROOT / "generated" / "p2635_s1585_toe_neural_universe_empirical_signature_audit.json"
MD = ROOT / "generated" / "p2635_s1585_toe_neural_universe_empirical_signature_audit.md"


class P2635ToeNeuralUniverseEmpiricalSignatureAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_is_present_and_nonempty(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        self.assertTrue(all(data["count"] > 0 for data in audit["patterns"].values()))
        self.assertIn("neural_universe_self_learning_content", audit["patterns"])
        self.assertIn("modern_physics_empirical_content", audit["patterns"])

    def test_toe_symptom_axes_identify_strongest_signs_and_blockers(self) -> None:
        axes = self.payload["toe_signature_axes"]
        self.assertEqual(len(axes), 8)
        score = self.payload["weighted_toe_signature_score"]
        self.assertGreater(score["toe_likeness_score_0_to_1_not_probability"], 0.4)
        self.assertLess(score["toe_likeness_score_0_to_1_not_probability"], 0.8)
        self.assertIn("single-kernel cross-sector unification signature", score["top_positive_axes"])
        self.assertIn("selector/source closure signature", score["top_blocking_axes"])

    def test_neural_self_learning_is_classified_as_conditional_not_proven(self) -> None:
        axes = {row["axis"]: row for row in self.payload["toe_signature_axes"]}
        learning = axes["self-learning / variational stationarity signature"]
        self.assertEqual(learning["closure_status"], "conditional_learning_proxy_source_not_derived")
        self.assertTrue(learning["computable_indicators"]["p2506_selector_postulated_not_derived"])
        neural = axes["neural architecture signature: positional encoding plus heavy-tailed attention"]
        self.assertEqual(neural["closure_status"], "architectural_analogy_positive_but_beta_criticality_unclosed")

    def test_modern_physics_test_matrix_is_falsifiable(self) -> None:
        matrix = self.payload["modern_physics_test_matrix"]
        self.assertEqual(len(matrix), 5)
        routes = [row["test_route"] for row in matrix]
        self.assertTrue(any("CMB/LSS" in route for route in routes))
        self.assertTrue(any("GW/PTA" in route for route in routes))
        self.assertTrue(all(row["pass_condition"] and row["failure_condition"] for row in matrix))

    def test_acceptance_preserves_nonclosure_flags(self) -> None:
        acceptance = self.payload["source_acceptance"]
        self.assertFalse(acceptance["accepts_toe_closure_or_self_learning_universe_proof"])
        self.assertIn("selector_source_closure_obtained", acceptance["failed_gates"])
        self.assertIn("empirical_modern_physics_confirmation_obtained", acceptance["failed_gates"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))

    def test_recommendation_and_docs_updated(self) -> None:
        self.assertIn("phase-frequency/node certificate", self.payload["recommended_next_honest_step"])
        self.assertIn("CMB/LSS", self.payload["recommended_next_honest_step"])
        self.assertIn("P2635/S1585", MD.read_text(encoding="utf-8"))
        self.assertIn("P2635/S1585", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2635/S1585", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
