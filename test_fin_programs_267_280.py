#!/usr/bin/env python3
"""Regression tests for FIN research Programs 267--280."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import unittest

import numpy as np

import fin_programs_255_266 as previous
import fin_programs_267_280 as programs


ROOT = Path(__file__).resolve().parent


class Programs267To280Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.results = json.loads(programs.RESULTS_PATH.read_text(encoding="utf-8"))

    def test_strict_generator_remains_positive(self) -> None:
        a, _ = previous.strict_operator()
        self.assertLess(np.linalg.norm(a @ np.ones(12)), 1e-12)
        self.assertGreater(np.linalg.eigvalsh(a).min(), -1e-12)

    def test_p267_unique_visible_measure_and_instability(self) -> None:
        result = self.results["P267"]
        self.assertEqual(len(result["visible_atomic_poles"]), 3)
        self.assertGreater(result["pole_separation_minimum"], 0.3)
        low = result["noise_stability"]["1e-06"][
            "median_max_relative_pole_error"
        ]
        high = result["noise_stability"]["0.001"][
            "median_max_relative_pole_error"
        ]
        self.assertGreater(high, low)

    def test_p268_category_certificate(self) -> None:
        result = self.results["P268"]
        self.assertEqual(result["chains_audited"], 50_000)
        self.assertLess(result["maximum_composition_residual"], 1e-12)
        self.assertGreater(result["minimum_reduced_positive_eigenvalue"], 0.0)

    def test_p269_minimal_current_povm(self) -> None:
        result = self.results["P269"]
        self.assertEqual(result["minimum_outcomes_by_normalization_bound"], 6)
        self.assertEqual(result["constructed_outcomes"], 6)
        self.assertEqual(result["response_rank"], 5)
        self.assertLess(result["povm_completeness_residual"], 1e-12)
        self.assertGreater(result["minimum_effect_eigenvalue"], 0.0)
        self.assertLess(result["exact_coefficient_reconstruction_residual"], 1e-12)

    def test_p270_rank_defect(self) -> None:
        rows = self.results["P270"]["rows"]
        self.assertEqual([row["context_size"] for row in rows], [6, 3, 1])
        self.assertTrue(all(row["lifted_rank"] < row["strict_rank"] for row in rows))
        self.assertTrue(all(row["relative_best_fit_residual"] > 0.69 for row in rows))

    def test_p271_quantum_information_ledger(self) -> None:
        result = self.results["P271"]
        self.assertLess(result["instrument_completeness_residual"], 1e-12)
        self.assertGreaterEqual(result["minimum_information_loss_nats"], -1e-12)
        self.assertLess(result["maximum_chain_decomposition_residual"], 1e-12)

    def test_p272_positive_shift_not_bridge(self) -> None:
        result = self.results["P272"]
        self.assertLess(result["legacy_minimum_mean_zero_eigenvalue"], 0.0)
        self.assertGreater(result["completed_minimum_mean_zero_eigenvalue"], 0.0)
        self.assertGreater(result["projective_strict_fingerprint_defect"], 0.4)
        self.assertGreater(result["normalized_self_energy_defect_at_z_0_2"], 0.8)

    def test_p273_no_false_external_admission(self) -> None:
        result = self.results["P273"]
        self.assertEqual(result["production_p241_manifests_found"], 0)
        self.assertEqual(result["template_event_rows"], [0, 0])
        self.assertFalse(result["p242_execution_authorized"])

    def test_p274_has_interior_noise_bias_optimum(self) -> None:
        result = self.results["P274"]
        optimum = result["optimal_h_on_frozen_grid"]
        grid = [row["h"] for row in result["grid"]]
        self.assertGreater(optimum, min(grid))
        self.assertLess(optimum, max(grid))
        self.assertLess(result["optimal_p95_relative_error"], 0.02)

    def test_p275_frozen_null_exclusion(self) -> None:
        result = self.results["P275"]
        self.assertEqual(result["null_count"], 400)
        self.assertEqual(result["false_positive_count"], 0)
        self.assertGreater(result["minimum_fingerprint_distance"], result["tolerance"])
        self.assertGreater(result["minimum_positive_certified_radius"], 0.0)

    def test_p276_two_time_theorem_witnesses(self) -> None:
        result = self.results["P276"]
        self.assertLess(result["static"]["relative_log_generator_defect"], 1e-12)
        self.assertGreater(
            result["shape_drift"]["relative_log_generator_defect"], 0.05
        )
        self.assertGreater(
            result["mixture_memory_model"]["relative_log_generator_defect"], 0.01
        )
        self.assertLess(
            result["uncalibrated_scale_clock_degeneracy_residual"], 1e-12
        )

    def test_p277_interventions_identify_frozen_panel(self) -> None:
        result = self.results["P277"]
        self.assertEqual(result["holdout_winner"], result["true_law"])
        self.assertEqual(result["design_rank_without_intervention"], 1)
        self.assertEqual(result["design_rank_with_intervention"], 2)
        rows = {row["law"]: row for row in result["candidate_rows"]}
        self.assertLess(
            rows["driven_additive"]["holdout_mse"],
            rows["relaxation"]["holdout_mse"] / 1e5,
        )

    def test_p278_nonlinear_benchmark_is_task_specific(self) -> None:
        result = self.results["P278"]
        self.assertEqual(result["control_count"], 60)
        self.assertEqual(result["task_count"], 4)
        self.assertEqual(result["tasks_on_which_fin_exceeds_all_controls"], 1)
        self.assertEqual(result["fin_percentiles"]["narma10_r2"], 1.0)
        self.assertLess(result["fin_percentiles"]["parity3_accuracy"], 0.1)

    def test_p279_conditioned_thermodynamic_identity(self) -> None:
        result = self.results["P279"]
        self.assertLess(
            result["maximum_gibbs_free_energy_identity_residual"], 1e-12
        )
        self.assertLess(result["information_to_energy_conversion_residual"], 1e-12)
        self.assertGreater(result["landauer_one_bit_lower_bound"], 0.0)

    def test_p280_product_torsor_independence(self) -> None:
        result = self.results["P280"]
        self.assertLess(result["scale_fingerprint_invariance_residual"], 1e-12)
        self.assertLess(result["orientation_sign_pair_residual"], 1e-12)
        self.assertEqual(len(result["rows"]), 4)

    def test_summary_contains_all_fourteen_programs(self) -> None:
        with programs.SUMMARY_PATH.open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 14)
        self.assertEqual(rows[0]["program"], "P267")
        self.assertEqual(rows[-1]["program"], "P280")

    def test_raw_table_cardinalities(self) -> None:
        expectations = {
            programs.STABILITY_PATH: 120,
            programs.POVM_PATH: 6,
            programs.RG_PATH: 3,
            programs.RESERVOIR_PATH: 61,
        }
        for path, expected in expectations.items():
            with path.open("r", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), expected, path.name)

    def test_all_five_figures_exist(self) -> None:
        expected = {
            "p267_p269_measure_povm.png",
            "p270_p272_rg_completion.png",
            "p274_p275_tomography_false_positive.png",
            "p277_p278_adaptation_reservoir.png",
            "p279_p280_conversion_torsors.png",
        }
        self.assertEqual(
            {path.name for path in programs.FIGURE_DIR.glob("*.png")}, expected
        )


if __name__ == "__main__":
    unittest.main()
