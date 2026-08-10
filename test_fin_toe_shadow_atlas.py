#!/usr/bin/env python3
"""Regression tests for the counterfactual FIN ToE shadow atlas."""

from __future__ import annotations

import json
import unittest
from pathlib import Path

import numpy as np

from fin_toe_shadow_atlas_analysis import (
    exact_operator_shadow_audit,
    laplacian,
    radial_matrix,
    strict_profile,
)


ROOT = Path(__file__).resolve().parent


class ShadowAtlasTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = json.loads((ROOT / "FIN_TOE_Shadow_Atlas_Results.json").read_text(encoding="utf-8"))

    def test_common_spectrum_reconstruction(self) -> None:
        errors = self.data["operator_audit"]["common_spectrum_reconstruction_max_errors"]
        self.assertTrue(all(value < 1e-12 for value in errors.values()))

    def test_live_operator_audit_matches_serialized_result(self) -> None:
        live = exact_operator_shadow_audit()
        frozen = self.data["operator_audit"]
        self.assertAlmostEqual(live["strict_row_sum_s"], frozen["strict_row_sum_s"], places=14)
        self.assertEqual(live["distinct_eigenvalue_count"], frozen["distinct_eigenvalue_count"])
        self.assertTrue(
            np.allclose(live["eigenvalues"], frozen["eigenvalues"], atol=1e-14, rtol=0)
        )

    def test_strict_matrix_uses_zero_diagonal(self) -> None:
        w = radial_matrix(strict_profile)
        a = laplacian(w)
        self.assertTrue(np.array_equal(np.diag(w), np.zeros(12)))
        self.assertLess(float(np.ptp(w.sum(axis=1))), 1e-14)
        self.assertLess(float(np.max(np.abs(a.sum(axis=1)))), 1e-14)

    def test_unitary_inverse_is_on_nonaliased_branch(self) -> None:
        premise = self.data["operator_audit"]["unitary_reconstruction_premise"]
        self.assertTrue(premise["known_zero_mode"])
        self.assertLess(premise["t_lambda_max"], premise["principal_branch_bound_pi"])

    def test_functional_algebra_is_commutative(self) -> None:
        audit = self.data["operator_audit"]
        self.assertEqual(audit["functional_algebra_dimension"], 7)
        self.assertEqual(audit["full_matrix_algebra_dimension"], 144)
        self.assertEqual(audit["commutative_dimension_deficit"], 137)
        self.assertTrue(all(value < 1e-12 for value in audit["all_functional_commutator_norms"].values()))

    def test_orientation_no_go_inside_functional_calculus(self) -> None:
        for row in self.data["operator_audit"]["symmetry_commutator_norms"].values():
            self.assertLess(row["translation"], 1e-12)
            self.assertLess(row["reflection"], 1e-12)

    def test_markov_semigroup(self) -> None:
        row = self.data["operator_audit"]["markov_checks"]
        self.assertGreater(row["minimum_entry"], 0)
        self.assertLess(row["maximum_row_sum_error"], 1e-12)
        self.assertLess(row["detailed_balance_uniform_error"], 1e-12)

    def test_resonance_energy_equivalence(self) -> None:
        self.assertLess(self.data["operator_audit"]["highest_resonance_equals_lowest_A_order_error"], 1e-12)

    def test_null_control(self) -> None:
        row = self.data["null_ensemble"]
        self.assertEqual(row["sample_count"], 10000)
        self.assertEqual(row["analytic_unitary_heat_green_property_rate_in_declared_class"], 1.0)
        self.assertEqual(row["generic_markov_semigroup_pass_rate"], 1.0)

    def test_legacy_strict_nonbridge(self) -> None:
        row = self.data["legacy_strict_fit"]
        self.assertGreater(row["best_relative_residual"], 0.9)
        self.assertGreater(row["canonical_beta_0_01_relative_residual"], 0.9)

    def test_shadow_ledger_has_no_physical_confirmation(self) -> None:
        rows = self.data["theory_shadow_ledger"]
        self.assertEqual(len(rows), 18)
        self.assertFalse(any(row["level"] == "L5" for row in rows))
        self.assertTrue(any(row["relation"].startswith("E0") for row in rows))
        self.assertTrue(any(row["relation"].startswith("E1") for row in rows))


if __name__ == "__main__":
    unittest.main(verbosity=2)
