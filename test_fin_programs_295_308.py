#!/usr/bin/env python3
"""Regression and scientific-boundary tests for FIN Programs P295--P308."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import subprocess
import unittest


ROOT = Path(__file__).resolve().parent
RESULTS = json.loads(
    (ROOT / "FIN_Programs_295_308_Results.json").read_text(encoding="utf-8")
)


class Programs295To308Tests(unittest.TestCase):
    def test_summary_contains_all_fourteen_programs(self) -> None:
        with (ROOT / "FIN_Programs_295_308_Summary.csv").open(
            encoding="utf-8"
        ) as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(
            [row["program"] for row in rows],
            [f"P{index}" for index in range(295, 309)],
        )

    def test_p295_multi_probe_is_full_rank_but_not_uniformly_stable(self) -> None:
        result = RESULTS["P295"]
        self.assertEqual(result["single_probe_jacobian_rank"], 6)
        self.assertEqual(result["four_probe_jacobian_rank"], 15)
        summaries = result["probe_summaries"]
        single = next(
            row
            for row in summaries
            if row["probe_count"] == 1 and row["noise_sigma"] == 0.0001
        )
        multiple = next(
            row
            for row in summaries
            if row["probe_count"] == 4 and row["noise_sigma"] == 0.0001
        )
        self.assertLess(
            multiple["median_maximum_pole_error"],
            single["median_maximum_pole_error"],
        )
        self.assertGreater(single["median_maximum_pole_error"], 0.1)

    def test_p296_formal_source_compiles_and_exact_checks_pass(self) -> None:
        result = RESULTS["P296"]
        self.assertEqual(result["exact_spd_rational_checks"], 400)
        self.assertEqual(result["lean_exit_code"], 0, result["lean_stderr"])
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

    def test_p297_complex_unitary_mesh_is_exact_before_quantization(self) -> None:
        result = RESULTS["P297"]
        self.assertEqual(result["n_modes"], 72)
        self.assertLess(result["exact_unitary_residual"], 1e-12)
        exact = result["rows"][0]
        self.assertLess(exact["unitary_reconstruction_residual"], 1e-12)
        self.assertLess(exact["maximum_effect_operator_error"], 1e-12)
        self.assertGreater(
            result["rows"][-1]["maximum_effect_operator_error"], 1e-5
        )

    def test_p298_flat_band_and_bernstein_obstructions(self) -> None:
        result = RESULTS["P298"]
        self.assertGreater(result["bernstein_mode_order_inversion_count"], 0)
        for row in result["rows"]:
            self.assertFalse(row["exact_unitary_equivalence_possible"])
            self.assertGreater(
                row["forced_flat_band_multiplicity"],
                row["strict_maximum_eigenvalue_multiplicity"],
            )
            self.assertGreater(row["optimized_projective_residual"], 0.01)

    def test_p299_memory_witness_separates_supplied_models(self) -> None:
        result = RESULTS["P299"]
        self.assertLess(abs(result["exact_homogeneous_cmi_nats"]), 1e-12)
        self.assertGreater(result["exact_hidden_memory_cmi_nats"], 1e-4)
        self.assertGreater(result["best_one_step_markov_prediction_tv"], 0.01)
        self.assertGreaterEqual(result["classification_accuracy"], 0.95)

    def test_p300_legacy_path_cone_does_not_contain_strict(self) -> None:
        result = RESULTS["P300"]
        self.assertEqual(result["sign_mismatch_count"], 4)
        self.assertLess(
            result["legacy_exact_laplace_mixture_quadrature_residual"], 1e-10
        )
        self.assertLess(
            result["strict_envelope_second_derivative_at_half_inflection"], 0.0
        )
        self.assertGreater(
            result["best_positive_fixed_phase_relative_residual"], 0.4
        )

    def test_p301_and_p306_external_gates_are_not_faked(self) -> None:
        self.assertFalse(RESULTS["P301"]["accepted_external_bundle"])
        self.assertFalse(
            RESULTS["P301"]["validator_241_executed_on_external_bundle"]
        )
        self.assertFalse(RESULTS["P306"]["accepted_hardware_record"])
        self.assertFalse(
            RESULTS["P306"]["independent_clock_and_detector_calibration"]
        )

    def test_p302_joint_nuisance_design_is_full_rank(self) -> None:
        result = RESULTS["P302"]
        self.assertEqual(result["optimal_fisher_rank"], 5)
        self.assertGreater(result["optimal_tested_h"], 0.1)
        self.assertLess(result["optimal_tested_h"], 1.8)
        self.assertLess(result["empirical_gradient_rmse"], 0.02)
        self.assertGreater(result["empirical_95_interval_coverage"], 0.85)

    def test_p303_small_ball_exponent_and_rank(self) -> None:
        result = RESULTS["P303"]
        self.assertEqual(result["fingerprint_jacobian_rank"], 5)
        self.assertEqual(result["local_asymptotic_exponent"], 5)
        self.assertGreater(result["gram_determinant"], 0.0)
        self.assertGreater(result["importance_fitted_exponent"], 4.0)
        self.assertLess(result["importance_fitted_exponent"], 6.5)

    def test_p304_clock_warp_matches_grid_but_not_holdout(self) -> None:
        result = RESULTS["P304"]
        self.assertGreater(result["minimum_clock_rate"], 0.0)
        self.assertLess(result["maximum_observed_time_propagator_error"], 1e-12)
        self.assertGreater(result["maximum_holdout_vertex_tv"], 1e-4)

    def test_p305_interventions_identify_supplied_coordinate_law(self) -> None:
        result = RESULTS["P305"]
        self.assertLess(result["passive_design_rank"], 6)
        self.assertEqual(result["joint_design_rank"], 6)
        self.assertLess(result["holdout_coordinate_rmse"], 1e-3)
        recovered = result["recovered_library_coefficients"]
        self.assertAlmostEqual(recovered["u"], 0.55, delta=0.02)
        self.assertAlmostEqual(recovered["lambda"], -0.32, delta=0.02)
        self.assertAlmostEqual(recovered["lambda^2"], -0.18, delta=0.02)

    def test_p307_scale_orbit_blocks_unchanged_role_transfer(self) -> None:
        result = RESULTS["P307"]
        self.assertTrue(result["strict_projective_orbit_exact"])
        self.assertEqual(result["clock_scale_fisher_rank"], 1)
        verdicts = result["legacy_role_verdicts"]
        self.assertFalse(verdicts["alpha_em_inverse"]["coordinate_scale_invariant"])
        self.assertFalse(
            verdicts["gravity_beta_power_N20"]["coordinate_scale_invariant"]
        )
        self.assertFalse(verdicts["sin2_theta_w"]["transferable"])

    def test_p308_pointing_does_not_create_a_strict_selector(self) -> None:
        result = RESULTS["P308"]
        for row in result["torsor_rows"]:
            self.assertEqual(row["invariant_global_sections"], 0)
            self.assertEqual(row["pointed_equivariant_maps"], 1)
        self.assertTrue(
            all(
                not certificate["certificate"]
                for certificate in result[
                    "historical_role_certificates"
                ].values()
            )
        )

    def test_raw_table_cardinalities(self) -> None:
        expected = {
            "FIN_Programs_295_308_MultiProbe_Stieltjes.csv": 96,
            "FIN_Programs_295_308_Optical_Mesh.csv": 4,
            "FIN_Programs_295_308_Spectral_Completion.csv": 2,
            "FIN_Programs_295_308_Process_Memory.csv": 120,
            "FIN_Programs_295_308_Sequential_Detector.csv": 25,
            "FIN_Programs_295_308_Rare_Event.csv": 4,
            "FIN_Programs_295_308_Adaptive_Discovery.csv": 880,
            "FIN_Programs_295_308_Legacy_Role_Invariance.csv": 7,
        }
        for filename, expected_count in expected.items():
            with (ROOT / filename).open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), expected_count, filename)

    def test_json_contains_no_nonstandard_numbers(self) -> None:
        text = (ROOT / "FIN_Programs_295_308_Results.json").read_text(
            encoding="utf-8"
        )
        self.assertNotIn("NaN", text)
        self.assertNotIn("Infinity", text)
        parsed = json.loads(text)
        self.assertTrue(math.isfinite(parsed["P303"]["gram_determinant"]))

    def test_all_five_figures_exist(self) -> None:
        expected = {
            "p295_p298_inverse_completion.png",
            "p297_p299_mesh_memory.png",
            "p302_p303_detector_rare_event.png",
            "p304_p305_clock_adaptation.png",
            "p300_p307_kernel_role_invariance.png",
        }
        actual = {
            path.name
            for path in (ROOT / "FIN_Programs_295_308_Figures").glob("*.png")
        }
        self.assertEqual(actual, expected)


if __name__ == "__main__":
    unittest.main()
