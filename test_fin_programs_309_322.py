#!/usr/bin/env python3
"""Regression and scientific-boundary tests for FIN Programs P309--P322."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import subprocess
import unittest


ROOT = Path(__file__).resolve().parent
RESULTS = json.loads(
    (ROOT / "FIN_Programs_309_322_Results.json").read_text(encoding="utf-8")
)


class Programs309To322Tests(unittest.TestCase):
    def test_summary_contains_all_fourteen_programs(self) -> None:
        with (ROOT / "FIN_Programs_309_322_Summary.csv").open(
            encoding="utf-8"
        ) as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(
            [row["program"] for row in rows],
            [f"P{index}" for index in range(309, 323)],
        )

    def test_p309_collision_has_quartic_information_modulus(self) -> None:
        result = RESULTS["P309"]
        self.assertAlmostEqual(
            result["collision_kl_power_exponent"], 4.0, delta=0.02
        )
        design = result["e_optimal_rows"]
        self.assertGreater(
            design[-1]["minimum_profiled_pole_information"],
            1.0e5 * design[0]["minimum_profiled_pole_information"],
        )
        resolutions = result["quarter_error_resolution_by_probe_count"]
        self.assertGreater(resolutions["1"], resolutions["4"])
        self.assertGreater(resolutions["4"], resolutions["16"])

    def test_p310_exact_blocks_and_formal_core(self) -> None:
        result = RESULTS["P310"]
        self.assertEqual(result["exact_rational_block_checks"], 300)
        self.assertEqual(result["lean_exit_code"], 0, result["lean_stderr"])
        self.assertFalse(result["mathlib_available_in_repository_environment"])
        lean = ROOT / ".elan/toolchains/leanprover--lean4---v4.28.0/bin/lean"
        process = subprocess.run(
            [str(lean), result["lean_source"]],
            cwd=ROOT,
            capture_output=True,
            text=True,
            timeout=120,
            check=False,
        )
        self.assertEqual(process.returncode, 0, process.stderr)

    def test_p311_loss_completed_measurement_remains_a_povm(self) -> None:
        result = RESULTS["P311"]
        self.assertEqual(result["all_tested_response_ranks"], [5])
        self.assertLess(result["maximum_completeness_residual"], 1e-12)
        self.assertGreater(result["minimum_effect_eigenvalue"], -1e-12)
        zero = next(
            row
            for row in result["grouped_results"]
            if row["mean_mode_loss"] == 0.0
            and row["phase_parameter_sigma"] == 0.0
        )
        worst = next(
            row
            for row in result["grouped_results"]
            if row["mean_mode_loss"] == 0.08
            and row["phase_parameter_sigma"] == 0.001
        )
        self.assertLess(zero["p95_maximum_effect_error"], 1e-12)
        self.assertGreater(worst["p95_maximum_effect_error"], 0.05)

    def test_p312_parent_is_cross_supported_but_target_dependent(self) -> None:
        rows = RESULTS["P312"]["rows"]
        for row in rows:
            size = row["carrier_size"]
            self.assertLess(row["signed_legacy_minimum_eigenvalue"], -1.0)
            self.assertGreater(row["normalized_cross_support"], 0.9)
            self.assertEqual(row["parent_rank"], size - 1)
            self.assertLess(abs(row["parent_minimum_eigenvalue"]), 1e-10)
            self.assertEqual(row["legacy_compression_residual"], 0.0)
            self.assertEqual(row["strict_compression_residual"], 0.0)
        self.assertGreater(
            rows[-1]["c12_trained_cubic_spectral_rmse"],
            rows[0]["c12_trained_cubic_spectral_rmse"],
        )

    def test_p313_positive_legacy_scale_mixture_has_wrong_tail(self) -> None:
        result = RESULTS["P313"]
        self.assertAlmostEqual(
            result["fitted_tail_log_slope"], -1.0, delta=0.01
        )
        self.assertAlmostEqual(
            result["strict_tail_log_slope"], -1.8, delta=0.01
        )
        self.assertGreater(
            result["positive_nnls_train_relative_residual"], 0.1
        )

    def test_p314_minimal_phase_extension_is_infinite_order(self) -> None:
        result = RESULTS["P314"]
        self.assertEqual(result["legacy_generated_phase_group_order"], 24)
        self.assertEqual(result["strict_common_angular_unit_radians"], "1/4000")
        self.assertEqual(result["strict_omega_integer_units"], 743)
        self.assertEqual(result["strict_phi_integer_units"], 650)
        self.assertEqual(result["integer_unit_gcd"], 1)

    def test_p315_operational_equivalence_fails_after_scale_fit(self) -> None:
        result = RESULTS["P315"]["results"]
        self.assertGreater(
            result["heat"]["finite_menu_operational_distance"], 0.2
        )
        self.assertGreater(
            result["wave"]["finite_menu_operational_distance"], 0.4
        )
        self.assertLess(
            result["wave"]["conservative_95pct_shot_bound"],
            result["heat"]["conservative_95pct_shot_bound"],
        )

    def test_p316_signed_resource_has_rigorous_nonzero_lower_bound(self) -> None:
        result = RESULTS["P316"]
        self.assertLess(result["fourth_forward_difference_at_zero"], 0.0)
        self.assertGreater(result["analytic_negative_mass_lower_bound"], 0.07)
        for row in result["grid_results"]:
            self.assertGreater(
                row["minimum_grid_negative_mass"],
                result["analytic_negative_mass_lower_bound"],
            )
            self.assertLess(row["moment_residual"], 1e-10)
        self.assertLess(
            result["grid_results"][-1]["minimum_grid_negative_mass"],
            result["grid_results"][0]["minimum_grid_negative_mass"],
        )

    def test_p317_and_p321_external_gates_remain_closed(self) -> None:
        for program in ("P317", "P321"):
            self.assertFalse(RESULTS[program]["admitted"])
            self.assertFalse(
                RESULTS[program]["one_shot_pipeline_authorized"]
            )
        self.assertFalse(RESULTS["P321"]["calibrated_clock_record"])

    def test_p318_tuned_dynamic_filter_beats_matched_static_filter(self) -> None:
        result = RESULTS["P318"]
        self.assertGreater(result["rmse_improvement_factor"], 1.2)
        self.assertTrue(result["best_dynamic_design"]["dynamic_filter"])
        self.assertLess(
            result["best_dynamic_design"]["median_gradient_rmse"], 0.007
        )

    def test_p319_exhausts_all_mode_order_chambers(self) -> None:
        result = RESULTS["P319"]
        self.assertEqual(result["possible_mode_orders"], 720)
        self.assertEqual(result["full_dimensional_feasible_chambers"], 720)
        self.assertEqual(result["full_rank_chamber_centers"], 720)
        self.assertGreater(result["minimum_center_gram_determinant"], 0.0)
        self.assertAlmostEqual(
            result["bonferroni_upper_bounds"]["720"],
            720 * result["p289_single_target_reference_probability"],
        )

    def test_p320_external_clock_slope_removes_one_rank_defect(self) -> None:
        result = RESULTS["P320"]
        self.assertEqual(
            result["uncalibrated_four_parameter_sensitivity_rank"], 3
        )
        self.assertEqual(result["uncalibrated_parameter_count"], 4)
        self.assertGreater(result["e_optimal_improvement_factor"], 2.0)
        self.assertEqual(len(result["best_design"]["times"]), 4)

    def test_p322_electroweak_role_certificate_fails_honestly(self) -> None:
        result = RESULTS["P322"]
        self.assertFalse(
            result["role_transfer_certificate"]["certificate_pass"]
        )
        self.assertTrue(
            all(
                not value
                for value in result["role_transfer_certificate"].values()
            )
        )
        self.assertAlmostEqual(
            result["legacy_benchmark"], 4.0 * math.log(2.0) / 12.0
        )
        self.assertGreater(
            result["candidate_stability"]["alpha_geo_over_carrier_size"][
                "range_across_carriers"
            ],
            0.1,
        )

    def test_raw_table_cardinalities(self) -> None:
        expected = {
            "FIN_Programs_309_322_Minimax_Stieltjes.csv": 114,
            "FIN_Programs_309_322_Lossy_Mesh.csv": 144,
            "FIN_Programs_309_322_Common_Parent.csv": 3,
            "FIN_Programs_309_322_Regular_Variation.csv": 183,
            "FIN_Programs_309_322_Operational_Distance.csv": 84,
            "FIN_Programs_309_322_Signed_Path.csv": 3,
            "FIN_Programs_309_322_Detector_Drift.csv": 560,
            "FIN_Programs_309_322_Spectral_Chambers.csv": 720,
            "FIN_Programs_309_322_Clock_Design.csv": 495,
            "FIN_Programs_309_322_Role_Naturality.csv": 15,
        }
        for filename, expected_count in expected.items():
            with (ROOT / filename).open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), expected_count, filename)

    def test_json_contains_no_nonstandard_numbers(self) -> None:
        text = (ROOT / "FIN_Programs_309_322_Results.json").read_text(
            encoding="utf-8"
        )
        self.assertNotIn("NaN", text)
        self.assertNotIn("Infinity", text)
        parsed = json.loads(text)
        self.assertTrue(
            math.isfinite(
                parsed["P320"]["best_design"][
                    "minimum_fisher_eigenvalue"
                ]
            )
        )

    def test_all_six_figures_exist(self) -> None:
        expected = {
            "p309_minimax_stieltjes.png",
            "p311_lossy_mesh.png",
            "p312_p313_parent_tail.png",
            "p315_p316_distance_negativity.png",
            "p318_p319_drift_chambers.png",
            "p320_p322_clock_role.png",
        }
        actual = {
            path.name
            for path in (ROOT / "FIN_Programs_309_322_Figures").glob("*.png")
        }
        self.assertEqual(actual, expected)


if __name__ == "__main__":
    unittest.main()
