#!/usr/bin/env python3
"""Regression and acceptance tests for FIN ST01--ST15."""

from __future__ import annotations

import hashlib
import json
import unittest
from pathlib import Path

import numpy as np

from fin_st01_st15_research import (
    CERTIFICATE,
    PREREGISTRATION,
    strict_operator,
    st02_common_spectrum_consistency,
)


ROOT = Path(__file__).resolve().parent


class ST01ST15Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = json.loads((ROOT / "FIN_ST01_ST15_Results.json").read_text(encoding="utf-8"))

    def test_all_fifteen_programs_present(self) -> None:
        self.assertTrue(all(f"ST{i:02d}" in self.data for i in range(1, 16)))

    def test_st01_composition_identity_and_bound(self) -> None:
        row = self.data["ST01"]
        self.assertLess(row["composition_identity_residual"], 1e-12)
        self.assertLessEqual(row["composition_norm_lhs"], row["composition_norm_bound_rhs"])
        self.assertTrue(row["negative_control_missing_instrument_rejected"])

    def test_st02_live_round_trip_and_negative_controls(self) -> None:
        _, a, _ = strict_operator()
        live = st02_common_spectrum_consistency(a)
        self.assertTrue(live["common_case_accepts"])
        self.assertTrue(live["altered_channel_rejected"])
        self.assertGreater(live["uncontrolled_alias_relative_error"], 0.5)
        self.assertTrue(all(error < 1e-12 for error in live["relative_reconstruction_errors"].values()))

    def test_st03_isospectral_and_null_controls(self) -> None:
        row = self.data["ST03"]
        self.assertLess(row["ensemble_summary"]["isospectral_rotations"]["maximum_isospectral_error"], 1e-12)
        self.assertGreater(row["ensemble_summary"]["isospectral_rotations"]["median_dihedral_residual"], 0.1)
        self.assertEqual(row["ensemble_summary"]["positive_radial"]["markov_m_matrix_rate"], 1.0)
        self.assertLess(row["ensemble_summary"]["signed_radial_kernel"]["psd_rate"], 0.5)

    def test_st04_algebra_completion_dimensions(self) -> None:
        dims = self.data["ST04"]["generated_algebra_dimensions"]
        self.assertEqual(dims["Cstar_A"], 7)
        self.assertEqual(dims["A_plus_distinct_vertex_labels"], 144)
        self.assertEqual(dims["A_plus_unoriented_anchor_distance"], 74)
        self.assertEqual(dims["A_plus_circulant_chiral_generator"], 12)

    def test_st05_continuum_nonuniqueness(self) -> None:
        row = self.data["ST05"]
        self.assertEqual(row["endpoint_reconstruction_error_N12"], 0.0)
        self.assertEqual(row["first_negative_weight_in_direct_integer_extension"], 8)
        self.assertLess(row["local_gap_scaling_exponent"], -1.8)
        self.assertGreater(row["dense_scaled_gap_scaling_exponent"], -0.3)

    def test_st06_strict_has_first_order_far_leakage(self) -> None:
        row = self.data["ST06"]
        strict = row["strict_rows"][0]["unitary_probability_outside_radius"]
        local = row["nearest_neighbor_control_rows"][0]["unitary_probability_outside_radius"]
        self.assertGreater(strict, 1e5 * local)
        self.assertGreater(row["minimum_direct_far_coupling"], 0)

    def test_st07_thermodynamic_identity_and_scale_orbit(self) -> None:
        row = self.data["ST07"]
        self.assertLess(row["gibbs_scale_orbit_error"], 1e-12)
        self.assertLess(row["relative_entropy_free_energy_identity_error"], 1e-12)
        self.assertEqual(len(row["axiom_removal_ledger"]), 5)

    def test_st08_gauge_receiver_identities(self) -> None:
        row = self.data["ST08"]
        self.assertLess(row["gauge_covariance_residual"], 1e-12)
        self.assertLess(row["holonomy_invariance_residual"], 1e-12)
        self.assertLess(row["continuity_equation_residual"], 1e-12)

    def test_st09_saturating_completion_is_coercive_in_sample(self) -> None:
        row = self.data["ST09"]
        self.assertGreater(row["sample_minimum_energy"], row["global_lower_bound"])
        self.assertGreater(row["large_amplitude_energy"], 100)
        self.assertIn("rho**2", row["small_density_series_through_rho3"])

    def test_st10_preserves_unresolved_exact_krein_boundary(self) -> None:
        row = self.data["ST10"]
        self.assertEqual(row["maximum_hamiltonian_structure_residual"], 0.0)
        self.assertFalse(row["positive_energy_sufficient_condition_ever_holds"])
        self.assertIsNotNone(row["sampled_transition_bracket"])
        self.assertIn("no_exact_Krein_closure", row["status"])

    def test_st11_exact_schur_and_truncated_scheme_dependence(self) -> None:
        row = self.data["ST11"]
        self.assertLess(row["exact_scheme_difference"], 1e-12)
        self.assertGreater(row["truncated_scheme_relative_difference"], 1e-4)

    def test_st12_blocks_role_transfer(self) -> None:
        row = self.data["ST12"]
        self.assertFalse(row["completion_endpoint_available"])
        self.assertFalse(row["role_invariance_theorem_available"])
        self.assertTrue(all("BLOCKED" in role["status"] or "REFUTED" in role["status"] or "NO_STRICT_TRANSFER" in role["status"] for role in row["roles"]))

    def test_st13_frozen_cross_channel_pipeline(self) -> None:
        row = self.data["ST13"]
        self.assertTrue(row["all_synthetic_same_A_channels_pass"])
        self.assertTrue(row["altered_generator_holdout_rejected"])
        self.assertTrue(all(error < 1e-12 for error in row["frozen_prediction_errors"].values()))

    def test_st14_preregistration_hash_and_holdout(self) -> None:
        prereg = json.loads(PREREGISTRATION.read_text(encoding="utf-8"))
        canonical = json.dumps(prereg["specification"], sort_keys=True, separators=(",", ":")).encode("utf-8")
        self.assertEqual(hashlib.sha256(canonical).hexdigest(), prereg["sha256"])
        row = self.data["ST14"]
        self.assertEqual(row["preregistration_sha256"], prereg["sha256"])
        self.assertEqual(row["holdout_false_accept_count"], 0)
        self.assertIn("failed_physical_prediction_requirement", row["status"])

    def test_st15_formal_certificate(self) -> None:
        certificate = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
        self.assertTrue(certificate["all_checks_pass"])
        self.assertTrue(all(certificate["checks"].values()))
        self.assertIn(11, certificate["Aut_Z12_unit_orbit"])

    def test_global_epistemic_boundary(self) -> None:
        boundary = self.data["epistemic_boundary"]
        for forbidden in ["QW-2191", "legacy-to-strict", "laboratory", "Theory of Everything"]:
            self.assertIn(forbidden, boundary)


if __name__ == "__main__":
    unittest.main(verbosity=2)
