"""Regression and scientific-boundary tests for Programs P323--P336."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import subprocess
import unittest


ROOT = Path(__file__).resolve().parent
RESULTS = json.loads(
    (ROOT / "FIN_Programs_323_336_Results.json").read_text(encoding="utf-8")
)


class Programs323To336Tests(unittest.TestCase):
    def test_summary_contains_all_fourteen_programs(self) -> None:
        with (ROOT / "FIN_Programs_323_336_Summary.csv").open(
            encoding="utf-8"
        ) as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(
            [row["program"] for row in rows],
            [f"P{index}" for index in range(323, 337)],
        )

    def test_p323_primal_dual_certificate_closes_grid_gap(self) -> None:
        result = RESULTS["P323"]
        self.assertAlmostEqual(
            result["exact_continuum_minimum_negative_mass_numeric"],
            0.40670633469225836,
            places=13,
        )
        self.assertLess(result["primal_dual_gap"], 1e-70)
        self.assertLess(result["maximum_moment_residual"], 1e-70)
        self.assertLessEqual(
            result["dual_q_squared_extremal_range"][1], 1.0 + 1e-13
        )

    def test_p324_same_input_has_distinct_targets(self) -> None:
        result = RESULTS["P324"]
        self.assertEqual(result["same_legacy_input_target_count"], 3)
        self.assertTrue(result["targets_pairwise_distinct"])
        self.assertGreater(
            result["maximum_target_distance_from_frozen"], 1.0
        )

    def test_p325_phase_unit_is_refreeze_lattice_not_internal_source(self) -> None:
        result = RESULTS["P325"]
        self.assertEqual(result["joint_additive_lattice_unit"], "1/4000")
        self.assertEqual(result["omega_coordinate"], 743)
        self.assertEqual(result["phi_coordinate"], 650)
        self.assertIn("[Refuted]", result["status"])

    def test_p326_lossy_protocol_remains_separating(self) -> None:
        best = RESULTS["P326"]["best_protocols"]
        self.assertGreater(best["heat"]["robust_tv_5pct"], 0.15)
        self.assertGreater(best["wave"]["robust_tv_5pct"], 0.30)
        self.assertLess(best["wave"]["conservative_95pct_shot_bound"], 100)

    def test_p327_all_chambers_have_exact_positive_witnesses(self) -> None:
        result = RESULTS["P327"]
        self.assertEqual(result["exact_distinct_orders"], math.factorial(6))
        self.assertTrue(result["all_exact_identities_hold"])
        self.assertEqual(result["global_minimum_weight_exact"], "7/50")
        self.assertEqual(result["mode_matrix_determinant"], "-3456*sqrt(3)")

    def test_p328_one_time_channel_reaches_full_diamond_distance(self) -> None:
        result = RESULTS["P328"]
        self.assertAlmostEqual(
            result["maximum_half_diamond_distance"], 1.0, places=12
        )
        self.assertIsNotNone(
            result["first_sampled_perfect_discrimination_time"]
        )
        self.assertIn("[Blocked]", result["status"])

    def test_p329_calibration_not_design_removes_scale_clock_defect(self) -> None:
        result = RESULTS["P329"]
        self.assertEqual(result["calibrated_nominal_rank"], 4)
        self.assertEqual(result["uncalibrated_rank"], 4)
        self.assertGreater(result["improvement_factor"], 3.0)
        self.assertTrue(
            result["best_design"]["all_sampled_clocks_monotone"]
        )

    def test_p330_and_p336_external_gates_remain_closed(self) -> None:
        for program in ("P330", "P336"):
            self.assertFalse(RESULTS[program]["admitted"])
            self.assertFalse(RESULTS[program]["one_shot_pipeline_authorized"])
            self.assertIn("[Blocked by external evidence]", RESULTS[program]["status"])

    def test_p331_tolerance_curve_contains_pass_and_fail(self) -> None:
        with (ROOT / "FIN_Programs_323_336_Photonic_Transfer.csv").open(
            encoding="utf-8"
        ) as handle:
            rows = list(csv.DictReader(handle))
        flags = [row["passes_declared_tv_budget"] == "True" for row in rows]
        self.assertIn(True, flags)
        self.assertIn(False, flags)
        self.assertGreater(
            RESULTS["P331"]["strict_nonzero_coupling_dynamic_range"], 40.0
        )

    def test_p332_negativity_does_not_select_interpretation(self) -> None:
        result = RESULTS["P332"]
        self.assertIn("[Refuted]", result["status"])
        self.assertIn("exp(i*pi)", result["phase_pi_realization"])

    def test_p333_cycle_commutes_but_noncyclic_carriers_obstruct(self) -> None:
        result = RESULTS["P333"]
        self.assertEqual(result["carrier_count"], 4)
        self.assertEqual(result["noncommuting_carriers"], 3)

    def test_p334_angle_package_is_explicitly_conditional(self) -> None:
        result = RESULTS["P334"]
        self.assertAlmostEqual(
            result["conditional_sin2_theta"],
            4.0 * math.log(2.0) / 12.0,
            places=14,
        )
        self.assertEqual(
            result["minimum_angle_only_additions"]["discrete_pointing_axiom"],
            1,
        )
        self.assertIn("[Refuted]", result["status"])

    def test_p335_five_single_omission_witnesses(self) -> None:
        result = RESULTS["P335"]
        self.assertEqual(result["configuration_count"], 32)
        self.assertTrue(result["all_single_omission_witnesses_exist"])
        self.assertEqual(len(result["single_omission_witnesses"]), 5)

    def test_formal_core_compiles(self) -> None:
        lean = Path(RESULTS["formal_verification"]["lean_binary"])
        completed = subprocess.run(
            [str(lean), "FIN_Programs_323_336_Formal_Core.lean"],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)

    def test_raw_table_cardinalities(self) -> None:
        expected = {
            "FIN_Programs_323_336_Continuum_Moment.csv": 18,
            "FIN_Programs_323_336_Source_Independent_Parent.csv": 3,
            "FIN_Programs_323_336_Phase_Provenance.csv": 6,
            "FIN_Programs_323_336_Lossy_Protocol.csv": 576,
            "FIN_Programs_323_336_Exact_Chambers.csv": 720,
            "FIN_Programs_323_336_Diamond_Distance.csv": 120,
            "FIN_Programs_323_336_Clock_Minimax.csv": 4368,
            "FIN_Programs_323_336_Photonic_Transfer.csv": 8,
            "FIN_Programs_323_336_Signed_Interpretations.csv": 6,
            "FIN_Programs_323_336_Carrier_Naturality.csv": 4,
            "FIN_Programs_323_336_EW_Axiom_Package.csv": 4,
            "FIN_Programs_323_336_Resource_Independence.csv": 32,
        }
        for filename, count in expected.items():
            with (ROOT / filename).open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), count, filename)

    def test_all_seven_figures_exist(self) -> None:
        expected = {
            "p323_continuum_extremizer.png",
            "p324_p325_parent_phase.png",
            "p326_lossy_protocol.png",
            "p327_exact_chambers.png",
            "p328_p329_channel_clock.png",
            "p331_p333_transfer_carriers.png",
            "p334_p335_axioms_independence.png",
        }
        actual = {
            path.name
            for path in (ROOT / "FIN_Programs_323_336_Figures").glob("*.png")
        }
        self.assertEqual(actual, expected)

    def test_json_contains_no_nonstandard_numbers(self) -> None:
        text = (ROOT / "FIN_Programs_323_336_Results.json").read_text(
            encoding="utf-8"
        )
        self.assertNotIn("NaN", text)
        self.assertNotIn("Infinity", text)


if __name__ == "__main__":
    unittest.main()
