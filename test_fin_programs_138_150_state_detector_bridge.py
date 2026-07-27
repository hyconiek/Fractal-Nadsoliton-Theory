#!/usr/bin/env python3
"""Regression tests for FIN Programs 138--150."""

import json
import math
import unittest
from pathlib import Path

import numpy as np

import fin_programs_138_150_state_detector_bridge as research


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "FIN_Programs_138_150_State_Detector_Bridge_Results.json"
PROTOCOL = ROOT / "FIN_Programs_138_150_PreData_Physical_Protocol.json"


class Programs138To150Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = json.loads(RESULTS.read_text())
        cls.protocol = json.loads(PROTOCOL.read_text())

    def test_01_release(self):
        self.assertEqual(self.data["release"], "10.14")

    def test_02_program_count(self):
        self.assertEqual(sum(k.startswith("program") for k in self.data), 13)

    def test_03_kms_does_not_select(self):
        p = self.data["program138_modular_kms_state_selection"]
        self.assertFalse(p["uniform_trace_selected"])
        self.assertFalse(p["strict_modular_gap_source_found"])

    def test_04_modular_hilbert_trace(self):
        p = self.data["program138_modular_kms_state_selection"]
        self.assertAlmostEqual(p["hilbert_gibbs_family"]["eta_at_theta_zero"], 17 / 9)

    def test_05_modular_target_gap(self):
        p = self.data["program138_modular_kms_state_selection"]
        self.assertAlmostEqual(
            p["hilbert_gibbs_family"]["theta_at_eta_9_over_5"], math.log(2)
        )

    def test_06_entropy_sector_counting(self):
        p = self.data["program139_maximum_entropy_reference"]
        self.assertAlmostEqual(p["reference_choices"]["sector_counting"]["eta"], 9 / 5)

    def test_07_entropy_hilbert(self):
        p = self.data["program139_maximum_entropy_reference"]
        self.assertAlmostEqual(p["reference_choices"]["Hilbert_microstate"]["eta"], 17 / 9)

    def test_08_entropy_not_reference_free(self):
        self.assertFalse(
            self.data["program139_maximum_entropy_reference"][
                "reference_free_uniform_state_found"
            ]
        )

    def test_09_morita_not_invariant(self):
        self.assertFalse(
            self.data["program140_morita_stability"][
                "normalized_Hilbert_trace_Morita_invariant"
            ]
        )

    def test_10_morita_uniform_amplification(self):
        self.assertAlmostEqual(
            self.data["program140_morita_stability"][
                "uniform_sector_eta_recovered_by_amplification"
            ],
            9 / 5,
        )

    def test_11_interval_fft_formal(self):
        p = self.data["program141_validated_interval_fft"]
        self.assertTrue(p["formal_interval_arithmetic"])
        self.assertTrue(p["target_window_fully_covered"])

    def test_12_interval_fft_positive(self):
        p = self.data["program141_validated_interval_fft"]
        self.assertTrue(p["all_symbol_lower_bounds_positive"])
        self.assertTrue(p["all_cells_compatible_with_fractional_prediction"])

    def test_13_interval_fft_honest_tightness(self):
        p = self.data["program141_validated_interval_fft"]
        self.assertFalse(p["tight_2_3_percent_claim_formally_recovered"])
        self.assertGreater(p["maximum_relative_remainder_upper"], 1)

    def test_14_continued_fraction_prefix(self):
        self.assertEqual(
            self.data["program142_diophantine_discrepancy"][
                "certified_continued_fraction_terms"
            ],
            [0, 16, 1, 10, 2, 67, 2, 2, 5, 1, 2, 928],
        )

    def test_15_discrepancy_bound(self):
        p = self.data["program142_diophantine_discrepancy"]
        self.assertTrue(p["all_observed_rows_below_theorem_bound"])
        self.assertLess(p["maximum_observed_to_bound_ratio"], 1)

    def test_16_no_global_power_claim(self):
        self.assertFalse(
            self.data["program142_diophantine_discrepancy"][
                "all_scale_power_remainder_proved"
            ]
        )

    def test_17_sobolev_threshold(self):
        self.assertEqual(
            self.data["program143_weighted_wave_estimates"]["threshold"], 0.5
        )

    def test_18_no_wave_smoothing(self):
        self.assertFalse(
            self.data["program143_weighted_wave_estimates"]["wave_smoothing"]
        )

    def test_19_detector_smooth_limit(self):
        p = self.data["program144_detector_resolution"]
        self.assertTrue(p["smooth_resolution_independent_limit"])
        self.assertFalse(p["delta_resolution_independent_limit"])

    def test_20_detector_quadratic_check(self):
        slope = self.data["program144_detector_resolution"][
            "observed_small_sigma_error_slope"
        ]
        self.assertAlmostEqual(slope, 2, delta=0.02)

    def test_21_joint_rank(self):
        p = self.data["program145_joint_identifiability"]
        self.assertEqual(p["joint_rank"], 3)
        self.assertEqual(p["one_time_rank"], 1)

    def test_22_joint_rmse(self):
        p = self.data["program145_joint_identifiability"]["synthetic_RMSE"]
        self.assertLess(p["alpha"], 0.003)
        self.assertLess(p["epsilon"], 0.001)
        self.assertLess(p["rho"], 0.02)

    def test_23_iqr_target(self):
        self.assertAlmostEqual(
            self.data["program146_calibration_invariant_observables"][
                "target_value_for_alpha_4_over_5"
            ],
            5 / 4,
        )

    def test_24_iqr_simulation(self):
        p = self.data["program146_calibration_invariant_observables"]["simulation"]
        self.assertAlmostEqual(p["mean_without_detector_noise"], 1.25, delta=0.02)
        self.assertAlmostEqual(p["mean_with_relative_detector_noise"], 1.25, delta=0.02)

    def test_25_reflection_witness_odd(self):
        p = self.data["program147_preparation_resource_theory"]["finite_checks"]
        self.assertLess(p["RCR_plus_C_norm"], 1e-12)

    def test_26_resource_mixture_free(self):
        p = self.data["program147_preparation_resource_theory"]["finite_checks"]
        self.assertAlmostEqual(p["M_symmetric_mixture"], 0, places=12)

    def test_27_resource_branches(self):
        p = self.data["program147_preparation_resource_theory"]["finite_checks"]
        self.assertAlmostEqual(p["M_rho_plus"], 1, places=12)
        self.assertAlmostEqual(p["M_rho_minus"], 1, places=12)
        self.assertAlmostEqual(p["Lambda_plus"] + p["Lambda_minus"], 0, places=12)

    def test_28_resource_bound(self):
        self.assertTrue(
            self.data["program147_preparation_resource_theory"]["finite_checks"][
                "all_mixing_rows_respect_witness_bound"
            ]
        )

    def test_29_selector_open(self):
        self.assertFalse(
            self.data["program147_preparation_resource_theory"][
                "QW_2191_discharged"
            ]
        )

    def test_30_variational_sector_solution(self):
        p = self.data["program148_coupled_state_damping_action"]
        self.assertAlmostEqual(p["sector_counting_solution"]["eta"], 9 / 5)

    def test_31_variational_hilbert_solution(self):
        p = self.data["program148_coupled_state_damping_action"]
        self.assertAlmostEqual(p["Hilbert_microstate_solution"]["eta"], 17 / 9)

    def test_32_variational_source_open(self):
        self.assertFalse(
            self.data["program148_coupled_state_damping_action"][
                "target_free_uniform_selection"
            ]
        )

    def test_33_bridge_control_commutes(self):
        p = self.data["program149_completion_map_diagram"]
        self.assertTrue(p["commutes_when_oscillatory_numerator_is_shared"])
        self.assertLess(p["shared_phase_frequency_control_residual"], 1e-14)

    def test_34_actual_bridge_not_commuting(self):
        p = self.data["program149_completion_map_diagram"]
        self.assertFalse(p["commutes_for_actual_canonical_legacy_parameters"])
        self.assertGreater(
            p["actual_canonical_legacy_to_projective_strict_relative_residual"],
            0.4,
        )

    def test_35_no_role_transfer(self):
        p = self.data["program149_completion_map_diagram"]
        self.assertFalse(p["full_bridge_complete"])
        self.assertFalse(p["role_transfer_started"])

    def test_36_protocol_digest(self):
        self.assertEqual(
            research.canonical_digest(self.protocol["core"]),
            self.protocol["canonical_core_sha256"],
        )

    def test_37_protocol_result_digest(self):
        self.assertEqual(
            self.data["program150_predata_protocol"]["canonical_core_sha256"],
            self.protocol["canonical_core_sha256"],
        )

    def test_38_protocol_no_external_data(self):
        self.assertFalse(self.protocol["core"]["external_data_admitted"])
        self.assertFalse(
            self.data["program150_predata_protocol"]["external_data_admitted"]
        )

    def test_39_protocol_power(self):
        p = self.data["program150_predata_protocol"]["synthetic_power_audit"]
        self.assertLessEqual(p["synthetic_false_positive_rate"], 0.01)
        self.assertGreaterEqual(p["synthetic_power"], 0.99)

    def test_40_thirteen_figures(self):
        self.assertEqual(
            len(list((ROOT / "FIN_Programs_138_150_State_Detector_Bridge_Figures").glob("*.png"))),
            13,
        )

    def test_41_fractional_circulant_psd(self):
        eig = np.linalg.eigvalsh(research.fractional_circulant())
        self.assertGreaterEqual(eig.min(), -1e-12)

    def test_42_global_closures_not_claimed(self):
        blocked = self.data["global_verdict"]["closures_not_claimed"]
        self.assertIn("canonical uniform sector trace", blocked)
        self.assertIn("QW-2191 selector closure", blocked)
        self.assertIn("full legacy-to-strict completion", blocked)
        self.assertIn("external experimental validation", blocked)


if __name__ == "__main__":
    unittest.main(verbosity=2)
