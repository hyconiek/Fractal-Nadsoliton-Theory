#!/usr/bin/env python3
"""Regression and boundary tests for FIN Programs P281--P294."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import subprocess
import unittest


ROOT = Path(__file__).resolve().parent
RESULTS = json.loads(
    (ROOT / "FIN_Programs_281_294_Results.json").read_text(encoding="utf-8")
)


class Programs281To294Tests(unittest.TestCase):
    def test_summary_contains_all_fourteen_programs(self) -> None:
        with (ROOT / "FIN_Programs_281_294_Summary.csv").open(
            encoding="utf-8"
        ) as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual([row["program"] for row in rows],
                         [f"P{index}" for index in range(281, 295)])

    def test_p281_bounded_fit_prevents_runaway_not_ill_conditioning(self) -> None:
        summaries = RESULTS["P281"]["summaries"]
        self.assertLess(
            summaries["0.001"]["p95_maximum_relative_pole_error"], 1.0
        )
        self.assertGreater(
            summaries["1e-06"]["median_maximum_relative_weight_error"], 0.5
        )
        self.assertEqual(RESULTS["P281"]["frozen_model_order"], 3)

    def test_p282_lean_kernel_record(self) -> None:
        result = RESULTS["P282"]
        self.assertTrue(result["kernel_exit_zero"])
        self.assertIn("(13 : Rat)/5", result["stdout"])

    def test_p282_lean_source_still_compiles(self) -> None:
        lean = ROOT / RESULTS["P282"]["lean_runtime"]
        process = subprocess.run(
            [str(lean), "FIN_Programs_281_294_Schur_Core.lean"],
            cwd=ROOT,
            capture_output=True,
            text=True,
            timeout=120,
            check=False,
        )
        self.assertEqual(process.returncode, 0, process.stderr)

    def test_p283_nai_mark_and_detector_rank(self) -> None:
        result = RESULTS["P283"]
        self.assertLess(result["isometry_residual"], 1e-12)
        self.assertEqual(result["calibrated_response_rank"], 5)
        self.assertGreater(result["minimum_effect_eigenvalue"], 0.0)
        self.assertLess(
            result["confusion_column_stochastic_residual"], 1e-14
        )

    def test_p284_complement_restores_rank_not_shape(self) -> None:
        for row in RESULTS["P284"]["rows"]:
            self.assertEqual(row["completed_rank"], 11)
            self.assertLess(
                row["completed_best_fit_residual"],
                row["bare_best_fit_residual"],
            )
            self.assertGreater(
                row["completed_projective_fingerprint_distance"], 0.1
            )

    def test_p285_multitime_ledger(self) -> None:
        result = RESULTS["P285"]
        self.assertLess(result["maximum_decomposition_residual"], 1e-12)
        self.assertLess(result["maximum_chain_record_residual"], 1e-12)
        self.assertTrue(
            all(value >= -1e-12 for value in result["minimum_losses"].values())
        )

    def test_p286_single_atom_fails_before_bridge(self) -> None:
        result = RESULTS["P286"]
        self.assertFalse(result["positive_generator"])
        self.assertLess(result["minimum_mean_zero_eigenvalue"], -1.0)
        self.assertIsNone(result["strict_projective_fingerprint_defect"])

    def test_p287_external_gate_remains_closed(self) -> None:
        result = RESULTS["P287"]
        self.assertEqual(result["production_manifest_count"], 0)
        self.assertEqual(result["template_event_rows"], [0, 0])
        self.assertFalse(result["p242_authorized"])

    def test_p288_detector_optimum_is_interior(self) -> None:
        result = RESULTS["P288"]
        grid = result["grid"]
        optimum = result["optimal_tested_h"]
        self.assertGreater(optimum, grid[0]["h"])
        self.assertLess(optimum, grid[-1]["h"])
        self.assertLess(result["optimal_median_relative_error"], 0.02)

    def test_p289_importance_sampling_is_stable_but_conditional(self) -> None:
        result = RESULTS["P289"]
        self.assertGreater(
            result["estimated_false_positive_probability"], 0.0
        )
        self.assertGreater(result["event_effective_sample_size"], 200.0)
        self.assertLess(
            result["proposal_sensitivity_max_to_min_ratio"], 1.5
        )
        self.assertEqual(result["direct_event_count"], 0)

    def test_p290_robustness_curve(self) -> None:
        accuracy = RESULTS["P290"]["classification_accuracy_by_noise"]
        self.assertGreater(accuracy["0.2"], 0.9)
        self.assertGreaterEqual(accuracy["0.05"], accuracy["0.35"])
        self.assertLess(accuracy["0.35"], 1.0)

    def test_p291_intervention_has_nonzero_information(self) -> None:
        result = RESULTS["P291"]
        self.assertEqual(
            result["named_designs"]["zero"][
                "information_minimum_eigenvalue"
            ],
            0.0,
        )
        self.assertGreater(
            result["best_design"]["information_minimum_eigenvalue"], 1.0
        )

    def test_p292_replication_is_finite_and_not_universal(self) -> None:
        result = RESULTS["P292"]
        self.assertEqual(result["replications"], 24)
        self.assertTrue(math.isfinite(result["median_advantage_over_best_control"]))
        self.assertGreater(result["primary_win_rate"], 0.0)
        self.assertLess(result["primary_win_rate"], 1.0)
        self.assertLess(result["p05_advantage_over_best_control"], 0.0)

    def test_p293_erasure_ledger(self) -> None:
        result = RESULTS["P293"]
        self.assertGreaterEqual(
            result["total_work"] + 1e-12,
            result["generalized_landauer_bound"],
        )
        self.assertGreater(result["total_entropy_production"], 0.0)
        self.assertLess(result["energy_balance_residual"], 1e-12)

    def test_p294_product_torsor_counts(self) -> None:
        result = RESULTS["P294"]
        self.assertEqual(result["equivariant_sections_from_trivial_input"], 0)
        self.assertEqual(result["equivariant_maps_from_resource_torsor"], 6)
        self.assertEqual(result["pointed_equivariant_maps"], 1)

    def test_raw_table_cardinalities(self) -> None:
        expected = {
            "FIN_Programs_281_294_Stieltjes_Recovery.csv": 180,
            "FIN_Programs_281_294_Detector_POVM.csv": 7,
            "FIN_Programs_281_294_RG_Completion.csv": 2,
            "FIN_Programs_281_294_Detector_Flux.csv": 18,
            "FIN_Programs_281_294_Mechanism_Classification.csv": 64,
            "FIN_Programs_281_294_Intervention_Design.csv": 244,
            "FIN_Programs_281_294_Reservoir_Replication.csv": 312,
        }
        for filename, count in expected.items():
            with (ROOT / filename).open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), count, filename)

    def test_json_contains_no_nonstandard_numbers(self) -> None:
        text = (ROOT / "FIN_Programs_281_294_Results.json").read_text(
            encoding="utf-8"
        )
        self.assertNotIn("NaN", text)
        self.assertNotIn("Infinity", text)

    def test_all_five_figures_exist(self) -> None:
        expected = {
            "p281_p283_recovery_povm.png",
            "p284_p286_rg_bridge.png",
            "p288_p289_detector_false_positive.png",
            "p290_p291_mechanism_intervention.png",
            "p292_p294_replication_thermo_torsor.png",
        }
        actual = {
            path.name
            for path in (ROOT / "FIN_Programs_281_294_Figures").glob("*.png")
        }
        self.assertEqual(actual, expected)


if __name__ == "__main__":
    unittest.main()
