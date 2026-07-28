#!/usr/bin/env python3
"""Regression tests for FIN Programs 255--266.

Passing this suite certifies deterministic finite calculations only.  It does
not certify physical realizability, an external experiment, a selector, or a
legacy-to-strict completion theorem.
"""

from __future__ import annotations

import csv
import json
import unittest

import numpy as np

import fin_programs_255_266 as programs


class Programs255To266Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.results = json.loads(
            programs.RESULTS_PATH.read_text(encoding="utf-8")
        )

    def test_strict_generator_is_dirichlet_positive(self) -> None:
        a, w = programs.strict_operator()
        self.assertLess(np.max(np.abs(a.sum(axis=1))), 1e-12)
        self.assertGreaterEqual(np.linalg.eigvalsh(a).min(), -1e-12)
        self.assertTrue(np.all(w[np.triu_indices(programs.N, 1)] > 0))

    def test_p255_exhausts_all_nontrivial_contexts(self) -> None:
        p255 = self.results["P255"]
        self.assertEqual(p255["contexts_audited"], 2**programs.N - 2)
        self.assertGreater(p255["minimum_hidden_principal_eigenvalue"], 0.0)

    def test_p255_stieltjes_certificates(self) -> None:
        p255 = self.results["P255"]
        self.assertGreaterEqual(
            p255["minimum_signed_derivative_eigenvalue_orders_0_to_4"],
            -1e-11,
        )
        self.assertLess(
            p255["maximum_measure_reconstruction_residual"], 1e-11
        )
        self.assertLess(p255["maximum_measure_mass_residual"], 1e-11)

    def test_p256_minimal_realization(self) -> None:
        p256 = self.results["P256"]
        self.assertEqual(p256["controllable_observable_rank"], 5)
        self.assertEqual(p256["minimal_stieltjes_realization_dimension"], 5)
        self.assertEqual(p256["visible_self_energy_poles"], 3)
        self.assertLess(
            p256["decoupled_two_mode_augmentation_residual"], 1e-12
        )

    def test_p257_context_composition(self) -> None:
        p257 = self.results["P257"]
        self.assertEqual(p257["chains_audited"], 20_000)
        self.assertLess(p257["maximum_composition_residual"], 1e-12)
        self.assertEqual(p257["identity_residual"], 0.0)

    def test_p258_closed_chiral_derivative(self) -> None:
        p258 = self.results["P258"]
        self.assertLess(p258["maximum_derivative_residual"], 1e-9)
        self.assertLess(p258["maximum_covariance_residual"], 1e-12)
        self.assertLessEqual(p258["maximum_norm_to_bound_ratio"], 1.0)

    def test_p259_strict_decimation_has_no_dimension_cycle(self) -> None:
        cardinalities = [
            row["cardinality"] for row in self.results["P259"]["levels"]
        ]
        self.assertEqual(cardinalities, [12, 6, 3, 1])
        self.assertTrue(
            all(
                distance > 1e-3
                for distance in self.results["P259"][
                    "successive_descriptor_distances"
                ]
            )
        )

    def test_p260_information_ledger(self) -> None:
        p260 = self.results["P260"]
        self.assertLess(p260["telescoping_residual"], 1e-12)
        self.assertLess(
            p260["maximum_record_chain_rule_residual"], 1e-12
        )
        self.assertGreater(p260["minimum_step_environment_loss"], 0.0)
        self.assertGreater(p260["minimum_loss_over_1000_random_pairs"], 0.0)

    def test_p261_amplitude_only_map_fails_positive_class(self) -> None:
        p261 = self.results["P261"]
        self.assertLess(
            p261[
                "legacy_amplitude_absorbed_generator_minimum_eigenvalue"
            ],
            -1.0,
        )
        self.assertGreater(p261["legacy_positive_resolvent_ray_pole"], 0.0)
        self.assertGreater(
            min(p261["projective_self_energy_shape_defect"].values()), 0.5
        )

    def test_p262_projective_then_scale_protocol(self) -> None:
        p262 = self.results["P262"]
        self.assertEqual(p262["invalid_reconstructions"], 0)
        self.assertGreaterEqual(p262["joint_pass_rate"], 0.95)
        self.assertLess(
            p262["scale_orbit_transition_invariance_residual"], 1e-12
        )
        self.assertIn("external experiment not executed", p262["status"])

    def test_p263_vertex_povm_no_go(self) -> None:
        p263 = self.results["P263"]
        self.assertEqual(p263["vertex_population_residual"], 0.0)
        self.assertLess(p263["opposite_current_residual"], 1e-12)
        self.assertEqual(p263["five_current_observable_rank"], 5)
        self.assertGreater(
            max(abs(value) for value in p263["harmonic_currents_rho_plus_d1_to_d5"]),
            1e-3,
        )

    def test_p264_false_positive_atlas(self) -> None:
        rates = self.results["P264"]["ensemble_pass_rates"]
        self.assertEqual(rates["strict_target"]["all_five"], 1.0)
        self.assertEqual(rates["random_positive_dense"]["all_five"], 0.0)
        self.assertEqual(
            rates["random_positive_circulant"]["stieltjes_memory"], 1.0
        )
        self.assertEqual(
            rates["random_positive_circulant"]["context_composition"], 1.0
        )

    def test_p265_identifiability_quotient(self) -> None:
        p265 = self.results["P265"]
        self.assertEqual(
            p265["exact_uniform_generator_vs_clock_degeneracy_residual"],
            0.0,
        )
        self.assertGreater(
            p265["scenarios"]["nonuniform_generator_change"][
                "projective_fingerprint_drift"
            ],
            0.05,
        )
        self.assertGreater(
            p265["scenarios"]["apparatus_memory"][
                "two_step_semigroup_defect"
            ],
            0.05,
        )

    def test_p266_benchmark_does_not_claim_advantage(self) -> None:
        p266 = self.results["P266"]
        self.assertEqual(p266["control_count"], 120)
        self.assertEqual(
            p266["fin_percentile_among_controls"][
                "memory_capacity_delays_1_to_20"
            ],
            0.0,
        )
        self.assertIn("only if it outperforms", p266["verdict"])

    def test_summary_and_raw_table_cardinalities(self) -> None:
        with programs.SUMMARY_PATH.open(encoding="utf-8") as handle:
            summary = list(csv.DictReader(handle))
        with programs.FALSE_POSITIVE_PATH.open(encoding="utf-8") as handle:
            false_positive = list(csv.DictReader(handle))
        with programs.RESERVOIR_PATH.open(encoding="utf-8") as handle:
            reservoir = list(csv.DictReader(handle))
        self.assertEqual(len(summary), 12)
        self.assertEqual(len(false_positive), 401)
        self.assertEqual(len(reservoir), 121)


if __name__ == "__main__":
    unittest.main(verbosity=2)
