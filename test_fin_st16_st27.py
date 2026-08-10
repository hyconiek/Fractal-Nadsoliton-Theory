#!/usr/bin/env python3
"""Live acceptance tests for FIN Programs ST16--ST27."""

from __future__ import annotations

import hashlib
import json
import unittest
from pathlib import Path

import numpy as np

import fin_st16_st27_research as research


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "FIN_ST16_ST27_Results.json"


class ST16ST27Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.results = json.loads(RESULTS.read_text(encoding="utf-8"))
        cls.w, cls.a, _ = research.strict_operator()

    def test_all_twelve_programs_present(self) -> None:
        self.assertTrue(all(f"ST{i}" in self.results for i in range(16, 28)))

    def test_st16_jordan_reduction_and_transition(self) -> None:
        row = self.results["ST16"]
        self.assertLess(row["formula_vs_full_chain_maximum_residual"], 1e-12)
        self.assertLess(row["stationary_state_residual_inf"], 1e-10)
        self.assertGreater(row["selected_collision_speed"], 1.27)
        self.assertLess(row["selected_collision_speed"], 1.29)
        self.assertEqual(row["left_of_collision"]["unstable_count"], 1)
        self.assertEqual(row["right_of_collision"]["unstable_count"], 0)

    def test_st16_live_gauge_vector(self) -> None:
        state, _, poles, residues = research.final_memory_state(self.a)
        live = research.gauge_chain_data(self.a, state, poles, residues, 1.6)
        self.assertLess(live["gauge_residual"], 1e-12)
        self.assertLess(live["generalized_chain_residual"], 1e-12)

    def test_st17_invariant_algebra_no_go(self) -> None:
        row = self.results["ST17"]
        self.assertEqual(row["dihedral_invariant_real_symmetric_dimension"], 7)
        self.assertEqual(row["strict_functional_algebra_dimension"], 7)
        self.assertEqual(row["algebra_dimension_after_natural_candidates"], 7)
        self.assertEqual(row["full_matrix_algebra_dimension"], 144)

    def test_st18_live_positive_refinement_fiber(self) -> None:
        embedding = np.zeros((24, 12))
        for x in range(24):
            embedding[x, x % 12] = 1.0
        lift_a = research.lifted_laplacian_24(self.w, 0.17)
        lift_b = research.lifted_laplacian_24(self.w, 0.83)
        self.assertLess(np.linalg.norm(lift_a @ embedding - embedding @ self.a), 1e-12)
        self.assertLess(np.linalg.eigvalsh(lift_a)[0], 1e-12)
        self.assertGreater(np.linalg.norm(lift_a - lift_b), 1.0)

    def test_st19_preregistration_hash_and_separation(self) -> None:
        prereg = json.loads((ROOT / "FIN_ST19_Cross_Sector_Preregistration.json").read_text(encoding="utf-8"))
        payload = json.dumps(prereg["configuration"], sort_keys=True, separators=(",", ":")).encode("utf-8")
        self.assertEqual(hashlib.sha256(payload).hexdigest(), prereg["sha256"])
        row = self.results["ST19"]
        self.assertEqual(row["same_generator_false_rejections"], 0)
        self.assertEqual(row["altered_generator_false_accepts"], 0)
        self.assertLess(row["same_generator_maximum_score"], 1e-12)
        self.assertGreater(row["altered_generator_minimum_score"], 1e-4)

    def test_st20_interval_branches_and_projectors(self) -> None:
        row = self.results["ST20"]
        self.assertGreater(row["minimum_distinct_eigenvalue_gap_lower_bound"], 0.04)
        self.assertGreater(row["unitary_principal_branch_margin_lower_bound"], 2.0)
        self.assertGreater(row["wave_arccos_branch_margin_lower_bound"], 2.0)
        self.assertGreater(row["altered_rank_one_generator_relative_separation_lower_bound"], 0.02)

    def test_st21_bounds_dominate_computed_errors(self) -> None:
        row = self.results["ST21"]
        for item in row["finite_range_rows"]:
            self.assertLessEqual(item["actual_antipodal_amplitude"], item["series_remainder_bound"] + 1e-14)
        for item in row["long_range_rows_N192"]:
            self.assertLessEqual(item["actual_unitary_operator_error"], item["spectral_duhamel_bound"] + 1e-12)
            self.assertLessEqual(item["spectral_duhamel_bound"], item["absolute_row_tail_bound"] + 1e-12)
        self.assertEqual(row["first_negative_literal_extension_weight"], 8)

    def test_st22_operational_independence_countermodels(self) -> None:
        metrics = self.results["ST22"]["finite_countermodel_metrics"]
        self.assertGreater(metrics["preparation_choice_total_variation"], 0.9)
        self.assertGreater(metrics["measurement_choice_total_variation"], 0.9)
        self.assertEqual(metrics["same_marginals_different_composition_total_variation"], 0.5)
        self.assertLess(metrics["unit_scale_orbit_residual"], 1e-14)

    def test_st23_continuous_connection_no_go(self) -> None:
        row = self.results["ST23"]
        self.assertEqual(row["translation_invariant_antisymmetric_one_form_dimension"], 5)
        self.assertEqual(row["reflection_constraint_exact_rank"], 5)
        self.assertEqual(row["fully_dihedral_invariant_continuous_one_form_dimension"], 0)
        self.assertLess(row["strict_reflection_residual"], 1e-12)

    def test_st24_saturation_classification(self) -> None:
        rows = {item["asymptotic_growth_exponent_a"]: item for item in self.results["ST24"]["rows"]}
        self.assertEqual(rows[0.0]["predicted_verdict"], "coercive")
        self.assertEqual(rows[0.25]["predicted_verdict"], "coercive")
        self.assertEqual(rows[0.5]["predicted_verdict"], "unbounded_below")
        self.assertTrue(all(abs(item["quartic_density_coefficient"] + 0.5) < 1e-15 for item in rows.values()))

    def test_st25_does_not_claim_unreplayed_source(self) -> None:
        row = self.results["ST25"]
        if not row["replay_passed"]:
            self.assertIn("not_machine_checked", row["status"])
        self.assertEqual(hashlib.sha256((ROOT / row["source_file"]).read_bytes()).hexdigest(), row["source_sha256"])

    def test_st26_separates_first_and_second_order_dimensions(self) -> None:
        row = self.results["ST26"]
        self.assertIn("t/c", row["positive_scale_action"]["unitary_and_heat"])
        self.assertIn("sqrt(c)", row["positive_scale_action"]["wave_stiffness"])
        self.assertIn("T^-2", row["theorem"])

    def test_st27_adversarial_no_refit_protocol(self) -> None:
        row = self.results["ST27"]
        nonzero = [item for item in row["rows"] if item["rotation_size_epsilon"] > 0]
        self.assertLess(max(item["maximum_isospectral_error"] for item in nonzero), 1e-12)
        self.assertLess(max(item["maximum_sector_specific_refit_error"] for item in nonzero), 1e-14)
        self.assertEqual(row["first_sampled_epsilon_with_full_rejection"], 0.01)

    def test_global_epistemic_boundary(self) -> None:
        boundary = self.results["epistemic_boundary"]
        self.assertIn("do not discharge QW-2191", boundary)
        self.assertIn("laboratory data", boundary)
        self.assertIn("ToE", boundary)


if __name__ == "__main__":
    unittest.main(verbosity=2)
