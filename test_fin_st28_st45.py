#!/usr/bin/env python3
"""Independent local checks for FIN Release 10.57 (ST28--ST45)."""

from __future__ import annotations

import hashlib
import json
import unittest
from pathlib import Path

import numpy as np

import fin_st28_st45_research as research


ROOT = Path(__file__).resolve().parent
RESULTS = json.loads((ROOT / "FIN_ST28_ST45_Results.json").read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def preregistration_digest(path: Path) -> tuple[str, str]:
    packet = json.loads(path.read_text(encoding="utf-8"))
    canonical = json.dumps(packet["configuration"], sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(canonical).hexdigest(), packet["sha256"]


class TestST28ST45(unittest.TestCase):
    def test_01_program_inventory(self) -> None:
        self.assertEqual(RESULTS["programs"], "ST28-ST45")
        self.assertTrue(all(f"ST{index}" in RESULTS for index in range(28, 46)))

    def test_02_st28_interval_certificate(self) -> None:
        result = RESULTS["ST28"]
        lo, hi = result["collision_speed_interval"]
        self.assertLess(lo, hi)
        self.assertTrue(result["unique_positive_root_in_box"])
        self.assertGreater(result["root_derivative_lower_bound"], 0.5)
        self.assertTrue(all(row["strict_inclusion"] for row in result["radius_rows"]))

    def test_03_st29_live_bicommutant_dimensions(self) -> None:
        _, a, _ = research.strict_operator()
        uniform = np.eye(research.N) / research.N
        generic = np.diag(np.arange(1, research.N + 1, dtype=float))
        self.assertEqual(research.bicommutant_algebra_dimension([a, uniform]), (7, 22))
        self.assertEqual(research.bicommutant_algebra_dimension([a, generic]), (144, 1))
        rows = RESULTS["ST29"]["rows"]
        self.assertEqual([row["algebra_dimension_with_density_observable"] for row in rows], [7, 74, 144])
        self.assertEqual([row["orbit_size"] for row in rows], [1, 12, 24])

    def test_04_st30_associativity_is_nondiscriminating(self) -> None:
        result = RESULTS["ST30"]
        self.assertLess(max(row["composite_intertwining_residual"] for row in result["rows"]), 1e-12)
        self.assertGreater(result["maximum_pairwise_coherent_lift_difference"], 1.0)

    def test_05_st31_preregistration_and_negative_power_result(self) -> None:
        result = RESULTS["ST31"]
        computed, embedded = preregistration_digest(ROOT / "FIN_ST31_Finite_Count_Preregistration.json")
        self.assertEqual(computed, embedded)
        self.assertEqual(embedded, result["preregistration_sha256"])
        self.assertLess(result["common_false_rejection_rate"], 0.03)
        self.assertLess(result["altered_detection_power"], 0.10)
        self.assertIn("failed", result["status"])

    def test_06_st32_polynomial_bound(self) -> None:
        result = RESULTS["ST32"]
        self.assertGreater(result["analytic_convolution_constant"], result["sampled_convolution_constant_lower_estimate"])
        self.assertTrue(all(row["actual_strict_amplitude"] <= row["polynomial_bound"] for row in result["rows_t_0_1_C12"]))

    def test_07_st33_holonomy_classification(self) -> None:
        result = RESULTS["ST33"]
        self.assertAlmostEqual(result["reflection_fixed_fluxes_mod_2pi"][0], 0.0)
        self.assertAlmostEqual(result["reflection_fixed_fluxes_mod_2pi"][1], np.pi)
        self.assertLess(result["random_state_gradient_holonomy_residual_from_one"], 1e-12)

    def test_08_st34_nonzero_but_uniform_minimizer(self) -> None:
        result = RESULTS["ST34"]
        self.assertTrue(result["trial_state_negative_energy"])
        self.assertFalse(result["best_candidate"]["localized_by_IPR_threshold"])
        self.assertLess(result["best_candidate"]["density_uniformity_residual"], 1e-8)
        self.assertEqual(result["localized_stationary_candidate_count"], 0)

    def test_09_st35_no_machine_check_claim(self) -> None:
        result = RESULTS["ST35"]
        self.assertFalse(result["configured_toolchain"])
        self.assertFalse(result["network_install_attempted"])

    def test_10_st36_exact_scale_ranks(self) -> None:
        result = RESULTS["ST36"]
        self.assertEqual((result["known_time_rank"], result["known_time_nullity"]), (6, 1))
        self.assertEqual((result["unknown_time_rank"], result["unknown_time_nullity"]), (6, 4))

    def test_11_st37_preregistered_nuisance_test(self) -> None:
        result = RESULTS["ST37"]
        computed, embedded = preregistration_digest(ROOT / "FIN_ST37_Nuisance_Preregistration.json")
        self.assertEqual(computed, embedded)
        self.assertEqual(embedded, result["preregistration_sha256"])
        self.assertLess(result["common_holdout_false_rejection_rate"], 0.02)
        rows = {round(row["rotation_size"], 3): row for row in result["adversarial_rows"]}
        self.assertEqual(rows[0.01]["detection_power"], 1.0)

    def test_12_st38_joint_closure_fails(self) -> None:
        rows = RESULTS["ST38"]["rows"]
        self.assertEqual([row["density_algebra_dimension"] for row in rows], [7, 74, 74, 144])
        self.assertTrue(all(row["phase_gradient_holonomy_residual_from_one"] < 1e-12 for row in rows))
        self.assertGreater(rows[-1]["total_absolute_pair_current"], 0.9)

    def test_13_st39_relative_axiom_inventory(self) -> None:
        result = RESULTS["ST39"]
        self.assertEqual(result["axiom_count"], 9)
        self.assertEqual(set(result["axioms"]), set(result["removal_witnesses"]))

    def test_14_st40_live_operational_transport(self) -> None:
        _, a, _ = research.strict_operator()
        live = research.st40_carrier_neutral_information(a)
        self.assertLess(live["transported_record_residual"], 1e-12)
        self.assertGreater(live["fixed_instrument_total_variation"], 0.1)

    def test_15_st41_spectrum_is_not_complete(self) -> None:
        hierarchy = RESULTS["ST41"]["hierarchy"]
        spectral = next(row for row in hierarchy if row["level"] == "I1")
        operational = next(row for row in hierarchy if row["level"] == "I3")
        self.assertFalse(spectral["complete"])
        self.assertTrue(operational["complete"])

    def test_16_st42_realization_removal_obstruction(self) -> None:
        result = RESULTS["ST42"]
        self.assertEqual(len(result["removal_ledger"]), 5)
        self.assertIn("cannot be defined", result["theorem"])

    def test_17_st43_feedback_is_distinct_from_symmetry(self) -> None:
        result = RESULTS["ST43"]
        self.assertEqual(sum(result["sector_counts"]), 240)
        self.assertGreater(result["mean_selected_radius_mu_positive"], 0.5)
        self.assertLess(result["radius_mu_negative"], 1e-8)

    def test_18_st44_strict_degenerate_mode_is_partial(self) -> None:
        result = RESULTS["ST44"]
        self.assertEqual(result["eigenspace_dimension"], 2)
        self.assertEqual(result["selected_axis_density_algebra_dimension"], 74)
        self.assertEqual(result["selected_axis_orbit_size"], 12)

    def test_19_st45_channel_inequivalence(self) -> None:
        result = RESULTS["ST45"]
        self.assertEqual(result["common_zero_mode_dimension"], 1)
        self.assertLess(result["unitary_spectrum_unit_modulus_residual"], 1e-12)
        self.assertGreater(result["least_positive_mode_heat_contraction_from_one"], 0.3)

    def test_20_global_epistemic_boundary(self) -> None:
        boundary = RESULTS["epistemic_boundary"]
        for phrase in ("No result supplies laboratory data", "SI units", "ToE closure"):
            self.assertIn(phrase, boundary)


if __name__ == "__main__":
    unittest.main(verbosity=2)
