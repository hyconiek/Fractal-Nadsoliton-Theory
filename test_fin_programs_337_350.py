"""Regression and scientific-boundary tests for Programs P337--P350."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import subprocess
import unittest


ROOT = Path(__file__).resolve().parent
RESULTS = json.loads(
    (ROOT / "FIN_Programs_337_350_Results.json").read_text(encoding="utf-8")
)


class Programs337To350Tests(unittest.TestCase):
    def test_summary_contains_all_fourteen_programs(self) -> None:
        with (ROOT / "FIN_Programs_337_350_Summary.csv").open(
            encoding="utf-8"
        ) as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(
            [row["program"] for row in rows],
            [f"P{index}" for index in range(337, 351)],
        )

    def test_p337_exact_dual_lower_certificate_is_sharp(self) -> None:
        result = RESULTS["P337"]
        self.assertLess(
            float(result["maximum_bernstein_upper_bound_decimal"]), 1.0
        )
        self.assertGreater(result["certified_dual_lower_bound"], 0.40670632)
        self.assertLess(result["distance_from_certified_lower_bound"], 1e-8)
        self.assertLess(
            result["certified_objective_interval_width_from_moments"], 1e-20
        )

    def test_p338_full_oscillatory_resource_has_positive_bracket(self) -> None:
        result = RESULTS["P338"]
        self.assertEqual(result["negative_moment_orders"], [8, 9, 10, 11])
        self.assertLessEqual(result["continuum_feasible_dual_range"][1], 0.0)
        self.assertGreaterEqual(result["continuum_feasible_dual_range"][0], -1.0)
        self.assertGreater(result["numerical_bracket_width"], 0.0)
        self.assertLess(result["numerical_bracket_width"], 1e-3)
        self.assertGreater(
            result["continuum_dual_lower_bound"],
            result["attenuation_envelope_seven_moment_resource"] + 0.29,
        )
        self.assertLess(result["high_precision_primal_moment_residual"], 1e-80)

    def test_p339_refreeze_has_no_independent_holdout(self) -> None:
        result = RESULTS["P339"]
        self.assertFalse(result["independent_frozen_holdout_available"])
        self.assertTrue(result["continuous_best_candidate"]["all_pass"])
        self.assertGreater(result["score_improvement"], 0.03)
        self.assertGreater(result["local_pass_fraction"], 0.5)
        self.assertLess(result["local_pass_fraction"], 0.9)

    def test_p340_cycle_interpolant_fails_all_noncycle_carriers(self) -> None:
        result = RESULTS["P340"]
        self.assertEqual(result["cycle_interpolant_degree"], 6)
        self.assertLess(result["cycle_training_spectral_residual"], 1e-10)
        self.assertEqual(result["noncycle_transfer_failures"], 4)

    def test_p341_mesh_and_branch_certificate(self) -> None:
        result = RESULTS["P341"]
        self.assertEqual(result["mesh_two_mode_rotation_count"], 66)
        self.assertLess(result["ideal_reconstruction_residual"], 1e-12)
        self.assertAlmostEqual(
            result["single_time_branch_generator_distance"],
            2.0 * math.pi / result["witness_time"],
            places=12,
        )
        self.assertLess(result["single_time_branch_unitary_residual"], 1e-12)

    def test_external_gates_p342_p350_remain_closed(self) -> None:
        self.assertFalse(RESULTS["P342"]["admitted"])
        self.assertFalse(RESULTS["P350"]["admitted"])
        self.assertFalse(
            RESULTS["P342"]["continuous_refreeze_external_validation_authorized"]
        )
        self.assertFalse(RESULTS["P350"]["hardware_execution_performed"])

    def test_p343_parallel_use_ladder(self) -> None:
        earliest = RESULTS["P343"]["earliest_sampled_perfect_time_by_uses"]
        self.assertEqual(
            [earliest[str(index)] for index in range(1, 5)],
            [0.72, 0.36, 0.24, 0.18],
        )
        self.assertIn("[Blocked]", RESULTS["P343"]["status"])

    def test_p344_equispaced_minimizes_declared_gap(self) -> None:
        self.assertEqual(
            RESULTS["P344"]["best_maximum_gap_design"],
            "nonparametric_equispaced_design",
        )
        with (ROOT / "FIN_Programs_337_350_Nonparametric_Clock.csv").open(
            encoding="utf-8"
        ) as handle:
            rows = list(csv.DictReader(handle))
        separations = {
            row["design"]: float(row["two_clock_midpoint_separation"])
            for row in rows
        }
        self.assertLess(
            separations["nonparametric_equispaced_design"],
            separations["P329_quartic_Fisher_design"],
        )

    def test_p345_operator_invariants_are_carrier_dependent(self) -> None:
        coefficients = RESULTS["P345"][
            "tested_operator_invariant_coefficients_of_variation"
        ]
        self.assertTrue(all(value > 0.03 for value in coefficients.values()))
        self.assertEqual(RESULTS["P345"]["universal_formula_level_data"]["eta"], "9/5")

    def test_p346_markov_negativity_monotone(self) -> None:
        result = RESULTS["P346"]
        self.assertEqual(result["typed_resource_count"], 5)
        self.assertLessEqual(
            result["maximum_negative_mass_monotonicity_violation"], 1e-12
        )

    def test_p347_no_tested_coupling_law_is_accepted(self) -> None:
        result = RESULTS["P347"]
        self.assertEqual(result["candidate_count"], 4)
        self.assertEqual(result["accepted_internal_coupling_laws"], 0)
        self.assertGreater(result["minimum_phase_central_rmse"], 0.5)

    def test_p348_is_explicitly_conditional(self) -> None:
        result = RESULTS["P348"]
        self.assertEqual(result["one_parameter_jacobian_rank"], 1)
        self.assertEqual(result["observable_count"], 5)
        self.assertIn("[Refuted]", result["status"])
        self.assertAlmostEqual(
            result["conditional_angle"], 4.0 * math.log(2.0) / 12.0, places=14
        )

    def test_p349_conversion_frame_requires_rank_three(self) -> None:
        result = RESULTS["P349"]
        self.assertEqual(result["dimensionless_FIN_log_unit_jacobian_rank"], 0)
        self.assertEqual(result["best_conditioned_triple"]["rank"], 3)
        self.assertEqual(result["best_conditioned_triple"]["condition_number"], 1.0)
        self.assertGreater(result["full_rank_triples"], 0)

    def test_formal_core_compiles(self) -> None:
        verification = RESULTS["formal_verification"]
        lean = Path(verification["lean_binary"])
        completed = subprocess.run(
            [str(lean), "FIN_Programs_337_350_Formal_Core.lean"],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)

    def test_raw_table_cardinalities(self) -> None:
        expected = {
            "FIN_Programs_337_350_Interval_Certificate.csv": 18,
            "FIN_Programs_337_350_Refreeze_Sensitivity.csv": 258,
            "FIN_Programs_337_350_Natural_Completion.csv": 5,
            "FIN_Programs_337_350_Photonic_Compilation.csv": 11,
            "FIN_Programs_337_350_Parallel_Comb.csv": 388,
            "FIN_Programs_337_350_Nonparametric_Clock.csv": 2,
            "FIN_Programs_337_350_Carrier_Invariants.csv": 6,
            "FIN_Programs_337_350_Resource_Monotones.csv": 69,
            "FIN_Programs_337_350_Coupling_Laws.csv": 4,
            "FIN_Programs_337_350_EW_Falsification.csv": 5,
            "FIN_Programs_337_350_Metrology.csv": 20,
        }
        for filename, count in expected.items():
            with (ROOT / filename).open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), count, filename)

    def test_all_six_figures_exist(self) -> None:
        expected = {
            "p337_p338_moment_certificates.png",
            "p339_p340_refreeze_naturality.png",
            "p341_p343_photonic_comb.png",
            "p344_p345_clock_invariants.png",
            "p346_p347_resources_couplings.png",
            "p348_p349_conditional_metrology.png",
        }
        actual = {
            path.name
            for path in (ROOT / "FIN_Programs_337_350_Figures").glob("*.png")
        }
        self.assertEqual(actual, expected)

    def test_json_contains_no_nonstandard_numbers(self) -> None:
        text = (ROOT / "FIN_Programs_337_350_Results.json").read_text(
            encoding="utf-8"
        )
        self.assertNotIn("NaN", text)
        self.assertNotIn("Infinity", text)


if __name__ == "__main__":
    unittest.main()
