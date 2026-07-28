#!/usr/bin/env python3
"""Regression tests for the second-generation FIN puzzle atlas.

The tests certify only finite operator identities and deterministic numerical
experiments.  They do not certify a physical interpretation, a selector, or a
legacy-to-strict completion theorem.
"""

from __future__ import annotations

import itertools
import unittest

import fin_nadsoliton_second_generation_atlas as atlas


class SecondGenerationAtlasTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.certificates = atlas.numerical_certificates()

    def test_exact_combinatorial_coverage(self) -> None:
        identifiers = sorted(atlas.BASE)
        self.assertEqual(len(identifiers), 12)
        self.assertEqual(len(list(itertools.combinations(identifiers, 2))), 66)
        self.assertEqual(len(list(itertools.combinations(identifiers, 3))), 220)
        self.assertEqual(
            12 + 66 + 220 + len(atlas.LARGE_GROUPS),
            305,
        )

    def test_derived_object_identifiers_are_unique(self) -> None:
        identifiers = [row[0] for row in atlas.DERIVED_OBJECTS]
        self.assertEqual(len(identifiers), 15)
        self.assertEqual(len(set(identifiers)), 15)

    def test_generation_three_identifiers_are_unique(self) -> None:
        identifiers = [row[0] for row in atlas.GENERATION3]
        self.assertEqual(len(identifiers), 10)
        self.assertEqual(len(set(identifiers)), 10)

    def test_effective_action_stationarity(self) -> None:
        self.assertLess(
            self.certificates["effective_action_stationarity_residual"],
            1e-12,
        )

    def test_nested_schur_reduction_composes(self) -> None:
        self.assertLess(
            self.certificates["context_schur_associativity_residual"],
            1e-12,
        )

    def test_sampled_stieltjes_signed_derivatives_are_psd(self) -> None:
        minimum_eigenvalues = self.certificates[
            "operator_stieltjes_signed_derivative_min_eigenvalues_orders_0_to_4"
        ]
        self.assertEqual(set(minimum_eigenvalues), {"0", "1", "2", "3", "4"})
        self.assertGreaterEqual(min(minimum_eigenvalues.values()), -1e-12)

    def test_flux_current_identity(self) -> None:
        flux = self.certificates["flux_current_identity"]
        self.assertGreater(abs(flux["C_chi"]), 1e-3)
        self.assertLess(flux["residual"], 1e-9)

    def test_chiral_memory_is_nonzero_and_inversion_odd(self) -> None:
        chiral = self.certificates["chiral_memory_susceptibility"]
        self.assertGreater(chiral["operator_norm"], 1e-2)
        self.assertLess(chiral["inversion_odd_covariance_residual"], 1e-12)

    def test_calibration_orbit_invariance(self) -> None:
        self.assertLess(
            self.certificates["calibration_torsor_heat_invariance_residual"],
            1e-12,
        )

    def test_projective_spectrum_invariance(self) -> None:
        self.assertLess(
            self.certificates["projective_fingerprint_scaling_residual"],
            1e-12,
        )

    def test_information_contraction_audit(self) -> None:
        contraction = self.certificates["information_contraction_1000_pairs"]
        self.assertEqual(contraction["negative_count_below_minus_1e_12"], 0)
        self.assertGreaterEqual(contraction["minimum_loss"], -1e-12)

    def test_diffusion_diameter_contracts_on_grid(self) -> None:
        distances = self.certificates["maximum_diffusion_distance_squared"]
        values = [distances[key] for key in ["0.1", "0.25", "0.5", "1.0"]]
        self.assertTrue(all(left > right for left, right in zip(values, values[1:])))


if __name__ == "__main__":
    unittest.main(verbosity=2)
