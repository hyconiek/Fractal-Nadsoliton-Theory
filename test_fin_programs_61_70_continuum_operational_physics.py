#!/usr/bin/env python3
"""Regression tests for FIN Programs 61--70."""

import math
import unittest

import fin_programs_61_70_continuum_operational_physics as p


class Programs6170Tests(unittest.TestCase):
    def test_61_projective_defects_are_measured_not_assumed_zero(self):
        r = p.program61_continuum_functor()["families"]
        strict = r["strict_absolute"]
        legacy = r["legacy_absolute"]
        self.assertFalse(strict["operator_defect_monotone_decreasing"])
        self.assertTrue(legacy["green_defect_monotone_decreasing"])
        self.assertGreater(strict["rows"][-1]["green_projective_defect"], 0.2)
        self.assertGreater(legacy["rows"][-1]["green_projective_defect"], 0.2)

    def test_62_resistance_limit_and_green_divergence(self):
        r = p.program62_regularizer_independence()
        for name in ["strict", "legacy_absolute"]:
            last = r[name][-1]
            self.assertLess(last["resistance_relative_error"], 1e-6)
            self.assertAlmostEqual(last["delta_times_green_fro_norm"], 1.0, places=6)

    def test_63_schur_semigroup(self):
        r = p.program63_compression_semigroup()
        self.assertLess(
            r["direct_48_to_12_vs_sequential_48_to_24_to_12_relative_residual"],
            1e-12,
        )

    def test_64_krein_does_not_make_raw_diffusion_markov(self):
        r = p.program64_krein_signed_alternative()
        self.assertGreater(r["inertia"]["negative"], 0)
        self.assertLess(r["unitarity_residual"], 1e-12)
        self.assertGreater(r["raw_diffusion_operator_norm"], 1.0)
        self.assertLess(r["raw_diffusion_minimum_entry"], 0.0)
        self.assertGreater(r["indefinite_green_distance_negative_entries"], 0)

    def test_65_one_preparation_binary_discrimination(self):
        r = p.program65_operational_tomography()
        self.assertEqual(r["channel_difference_rank"], 11)
        self.assertGreater(r["one_localized_preparation_full_detector_TV"], 0.3)
        self.assertAlmostEqual(
            r["one_preparation_optimal_binary_difference"],
            r["one_localized_preparation_full_detector_TV"],
            places=12,
        )

    def test_66_causal_leading_orders(self):
        r = p.program66_causal_order()
        self.assertEqual(r["strict_all_distance_first_nonzero_power"], 1)
        self.assertEqual(r["strict_wave_leading_probability_exponent"], 2)
        self.assertEqual(r["nearest_neighbor_first_nonzero_power"], 6)
        self.assertEqual(r["nearest_neighbor_wave_leading_probability_exponent"], 12)

    def test_67_calibration_rank(self):
        r = p.program67_calibration_identifiability()
        self.assertEqual(r["rank_length_plus_time"], 2)
        self.assertEqual(r["rank_length_time_energy"], 3)
        self.assertLess(r["maximum_relative_recovery_error"], 1e-12)

    def test_68_chiral_formula_is_receiver_not_radial_source(self):
        r = p.program68_chiral_source_law()["results"]
        for row in r.values():
            self.assertEqual(row["radial_candidate"], 0.0)
            self.assertGreater(abs(row["directed_plus"]), 0.1)
            self.assertLess(row["sign_pair_residual"], 1e-12)
            self.assertLess(row["reflection_odd_residual"], 1e-12)

    def test_69_landauer_protocol(self):
        r = p.program69_landauer_protocol()
        for row in r["rows"]:
            self.assertLess(abs(row["first_law_residual"]), 1e-12)
            self.assertLess(abs(row["Clausius_equality_residual"]), 1e-12)
        self.assertAlmostEqual(
            r["limit_beta_gap_to_infinity"]["beta_work"], math.log(2.0), places=12
        )

    def test_70_blinded_heldout_challenge(self):
        r = p.program70_blinded_challenge()
        self.assertFalse(r["hidden_generator_in_candidate_set"])
        self.assertEqual(r["winner"], "strict")
        ranking = r["candidate_ranking"]
        self.assertLess(
            ranking[0]["heldout_mean_empirical_KL"],
            ranking[1]["heldout_mean_empirical_KL"],
        )


if __name__ == "__main__":
    unittest.main()
