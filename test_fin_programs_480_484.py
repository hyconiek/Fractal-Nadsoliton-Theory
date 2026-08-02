#!/usr/bin/env python3
"""Regression and adversarial tests for FIN local Programs P480 and P484."""

from __future__ import annotations

import ast
from copy import deepcopy
from fractions import Fraction
import json
from pathlib import Path
import unittest

import numpy as np

import fin_program_480_standalone_checker as checker


ROOT = Path(__file__).resolve().parent


class Program480StandaloneReplayTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.certificate = json.loads(
            (ROOT / "FIN_Program_480_Standalone_Certificate.json").read_text(
                encoding="utf-8"
            )
        )
        cls.saved_result = json.loads(
            (ROOT / "FIN_Program_480_Standalone_Check_Result.json").read_text(
                encoding="utf-8"
            )
        )

    def test_replay_recomputes_exact_acceptance(self) -> None:
        result = checker.replay(self.certificate)
        self.assertTrue(result["passed"])
        self.assertTrue(result["krawczyk_included"])
        self.assertTrue(result["normalizer_center_shifted_sylvester"])
        self.assertTrue(result["X3_center_shifted_sylvester"])
        self.assertGreater(Fraction(result["minimum_inclusion_margin"]), 0)
        self.assertLess(Fraction(result["maximum_contraction_row_sum"]), 1)
        self.assertEqual(result["exact_global_value_interval"], self.saved_result[
            "exact_global_value_interval"
        ])

    def test_zeroed_preconditioner_is_rejected(self) -> None:
        tampered = deepcopy(self.certificate)
        tampered["preconditioner"] = [
            ["0" for _ in row] for row in tampered["preconditioner"]
        ]
        result = checker.replay(tampered)
        self.assertFalse(result["passed"])
        self.assertFalse(result["krawczyk_included"])

    def test_checker_imports_only_the_declared_standard_library(self) -> None:
        source = (ROOT / "fin_program_480_standalone_checker.py").read_text(
            encoding="utf-8"
        )
        tree = ast.parse(source)
        imported: set[str] = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported.update(alias.name.split(".")[0] for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imported.add(node.module.split(".")[0])
        self.assertLessEqual(imported, {
            "__future__", "fractions", "json", "pathlib", "sys", "typing"
        })


class Program484ExactPhaseFaceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.result = json.loads(
            (ROOT / "FIN_Program_484_Results.json").read_text(encoding="utf-8")
        )
        cls.witness = np.load(ROOT / "FIN_Program_484_Phase_Face_Witness.npz")

    def test_exact_representation_checks_all_pass(self) -> None:
        self.assertTrue(all(cls_value for cls_value in self.result[
            "exact_representation_checks"
        ].values()))
        self.assertTrue(all(cls_value for cls_value in self.result[
            "exact_odd_dimension_counterexample"
        ]["checks"].values()))
        self.assertFalse(self.result["full_cone_optimizer_uniqueness_refuted"])
        self.assertFalse(self.result["exact_nonzero_optimal_segment_proved"])
        self.assertTrue(self.result["conditional_nonzero_optimal_segment_theorem"])
        self.assertEqual(self.result["numerical_phase_face_dimension_lower_bound"], 1)
        self.assertFalse(self.result["complete_optimal_face_dimension_proved"])

    def test_saved_phase_direction_is_causal_and_tangent(self) -> None:
        q_direction = self.witness["Q"]
        x3 = self.witness["X3"]
        delta = self.witness["Delta"]
        k_matrix = np.real(delta / 1j)
        self.assertLess(np.linalg.norm(q_direction + q_direction.T, np.inf), 2e-12)
        # Homogeneous three-comb causality for this P484 direction is the
        # block relation Q = Q0 direct-sum (-Q0), not a stochastic row-sum
        # constraint.
        self.assertLess(np.linalg.norm(q_direction[:4, 4:], np.inf), 2e-12)
        self.assertLess(np.linalg.norm(q_direction[4:, :4], np.inf), 2e-12)
        self.assertLess(
            np.linalg.norm(q_direction[:4, :4] + q_direction[4:, 4:], np.inf),
            2e-12,
        )
        residual = x3 @ q_direction @ x3 + k_matrix @ q_direction @ k_matrix / 4
        self.assertLess(np.linalg.norm(residual, np.inf), 2e-12)

    def test_phase_face_boundaries_are_preserved(self) -> None:
        text = self.result["boundary"]
        self.assertIn("does not prove nonuniqueness", text)
        self.assertIn("not classify the complete face", text)
        self.assertIn("supply physical units", text)
        metadata = self.result["metadata"]
        self.assertIn("QW-2191 remains open", metadata["selector_boundary"])
        self.assertIn("No legacy/strict substitution", metadata["kernel_boundary"])


class Program475ResourceBoundedEliminationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.result = json.loads(
            (ROOT / "FIN_Program_475_Results.json").read_text(encoding="utf-8")
        )

    def test_exact_trigonometric_field_reduction(self) -> None:
        field = self.result["exact_field_reduction"]
        self.assertTrue(field["exact"])
        self.assertEqual(field["field_degree"], 4)
        x = np.sin(np.pi / 8)
        self.assertLess(abs(8*x**4 - 8*x**2 + 1), 2e-15)
        self.assertLess(abs((1 - 2*x**2) - np.sin(np.pi/4)), 2e-15)
        self.assertLess(abs((3*x - 4*x**3) - np.sin(3*np.pi/8)), 2e-15)

    def test_affine_structure_and_pivot_are_recorded(self) -> None:
        inventory = self.result["system_inventory"]
        self.assertEqual(inventory["equations_before_field_relation"], 13)
        self.assertTrue(
            inventory["all_equations_affine_in_normalizer_coordinates_A_B_u"]
        )
        self.assertEqual(inventory["naive_cubic_Bezout_bound_over_Q_alpha"], 3**13)
        pivot = self.result["local_linear_pivot"]
        self.assertTrue(pivot["exact_pivot_determinant_is_nonzero_polynomial"])
        self.assertGreater(pivot["exact_pivot_determinant_term_count"], 0)

    def test_timeout_is_not_promoted_to_impossibility(self) -> None:
        self.assertFalse(self.result["minimal_polynomial_candidate_obtained"])
        self.assertIn("resource-bounded no-go", self.result["status"])
        self.assertIn("No claim is made that L lacks", self.result["not_proved"][0])
        attempt = self.result["bounded_lexicographic_attempt"]
        self.assertFalse(attempt["completed"])
        self.assertEqual(attempt["termination"], "wall_clock_timeout")
        self.assertEqual(attempt["wall_limit_seconds"], 45)


if __name__ == "__main__":
    unittest.main()
