#!/usr/bin/env python3
"""Regression tests for FIN Programs 151--164."""

from __future__ import annotations

import json
import math
import tempfile
import unittest
from pathlib import Path

import numpy as np

import fin_programs_151_164_axiomatic_operational_foundations as research


ROOT = Path(__file__).resolve().parent
RESULTS = json.loads(
    (ROOT / "FIN_Programs_151_164_Axiomatic_Operational_Results.json").read_text(
        encoding="utf-8"
    )
)
P = RESULTS["programs"]


class TestProgram151(unittest.TestCase):
    def test_doubled_fft(self):
        self.assertEqual(P["151"]["fft_length"], 32768)

    def test_window_cells(self):
        self.assertEqual(P["151"]["frequency_cells"], 100)

    def test_positive_and_compatible(self):
        self.assertTrue(P["151"]["all_symbol_lower_bounds_positive"])
        self.assertTrue(P["151"]["all_cells_fractional_compatible"])

    def test_no_false_tight_claim(self):
        self.assertFalse(P["151"]["sub_3_percent_certificate"])
        self.assertGreater(P["151"]["maximum_relative_remainder_upper"], 0.03)


class TestProgram152(unittest.TestCase):
    def test_cf_prefix(self):
        self.assertEqual(
            P["152"]["certified_terms"],
            [0, 16, 1, 10, 2, 67, 2, 2, 5, 1, 2, 928],
        )

    def test_kappa_positive_on_prefix(self):
        self.assertGreater(P["152"]["finite_prefix_kappa"]["nu_1"], 0)
        self.assertGreater(P["152"]["finite_prefix_kappa"]["nu_2"], 0)

    def test_global_not_claimed(self):
        self.assertFalse(P["152"]["all_scale_axiom_derived_for_theta"])
        self.assertTrue(P["152"]["finite_prefix_is_not_global_proof"])


class TestProgram153(unittest.TestCase):
    def test_orbits(self):
        self.assertEqual(P["153"]["structural_orbits"], [[2], [3], [5, 7, 11]])

    def test_simplex_dimension(self):
        self.assertEqual(P["153"]["simplex_dimension"], 2)

    def test_uniform_orbit_mass(self):
        np.testing.assert_allclose(
            P["153"]["uniform_sector_orbit_masses"], [0.2, 0.2, 0.6]
        )

    def test_no_unique_measure(self):
        self.assertFalse(
            P["153"]["normalization_additivity_naturality_select_unique_measure"]
        )


class TestProgram154(unittest.TestCase):
    def test_beta_exact_identity(self):
        self.assertAlmostEqual(P["154"]["beta_star"], math.log(2), places=14)
        self.assertAlmostEqual(P["154"]["beta_star"], RESULTS["constants"]["alpha_geo"] / 4)

    def test_uniform_central_state(self):
        np.testing.assert_allclose(
            P["154"]["central_weights_at_beta_star"], np.ones(5) / 5
        )

    def test_eta(self):
        self.assertAlmostEqual(P["154"]["eta_at_beta_star"], 9 / 5)

    def test_conditional_only(self):
        self.assertFalse(P["154"]["strict_core_derivation"])


class TestProgram155(unittest.TestCase):
    def test_cms_scale(self):
        rng = np.random.default_rng(9)
        x = research.stable_rvs(0.8, 10000, rng)
        self.assertTrue(np.all(np.isfinite(x)))

    def test_resource_monotone(self):
        self.assertEqual(P["155"]["monotone"], "M(rho_r)=|r|")

    def test_full_space_not_claimed(self):
        self.assertFalse(P["155"]["full_density_space_completeness"])

    def test_selector_not_sourced(self):
        self.assertFalse(P["155"]["strict_signed_preparation_source"])


class TestProgram156(unittest.TestCase):
    def test_asymptotic_sd_positive(self):
        self.assertGreater(P["156"]["predicted_zero_noise_slope_sd"], 0)

    def test_zero_noise_unbiased_in_mc(self):
        row = P["156"]["simulation"]["rows"][0]
        self.assertLess(abs(row["slope_bias"]), 0.01)

    def test_resolution_biases_down(self):
        rows = P["156"]["simulation"]["rows"]
        self.assertLess(rows[-1]["mean_slope"], rows[0]["mean_slope"])

    def test_deconvolution_warning_present(self):
        self.assertIn("unbounded", P["156"]["deconvolution_instability"])


class TestProgram157(unittest.TestCase):
    def test_constant_gain_identifiable(self):
        self.assertTrue(P["157"]["finite_tangent_test"][0]["alpha_identifiable"])

    def test_linear_gain_confounds(self):
        self.assertFalse(P["157"]["finite_tangent_test"][1]["alpha_identifiable"])

    def test_all_higher_polynomials_confound(self):
        self.assertTrue(
            all(
                not row["alpha_identifiable"]
                for row in P["157"]["finite_tangent_test"][1:]
            )
        )

    def test_rank_loss_exact(self):
        row = P["157"]["finite_tangent_test"][1]
        self.assertEqual(row["nuisance_rank"], row["augmented_rank"])


class TestProgram158(unittest.TestCase):
    def test_true_slope(self):
        self.assertAlmostEqual(P["158"]["true_T"], 1.25)

    def test_sd_decreases(self):
        sd = [row["empirical_sd_T"] for row in P["158"]["rows"]]
        self.assertTrue(all(a > b for a, b in zip(sd, sd[1:])))

    def test_coverage_reasonable(self):
        for row in P["158"]["rows"]:
            self.assertGreater(row["normal_95_coverage_T"], 0.90)
            self.assertLess(row["normal_95_coverage_T"], 0.99)

    def test_nonasymptotic_not_claimed(self):
        self.assertFalse(P["158"]["nonasymptotic_concentration_proved"])


class TestProgram159(unittest.TestCase):
    def test_protocol_was_frozen(self):
        self.assertTrue(P["159"]["protocol_rule_frozen_before_family_comparison"])

    def test_gaussian_not_rejected(self):
        row = P["159"]["rows"][0]
        self.assertLess(row["P150_fractional_decision_rate"], 0.05)

    def test_alpha_one_false_attribution(self):
        row = next(r for r in P["159"]["rows"] if r["family"] == "stable_alpha_1.0")
        self.assertGreater(row["P150_fractional_decision_rate"], 0.95)

    def test_apparatus_confounder_false_attribution(self):
        row = next(r for r in P["159"]["rows"] if r["family"] == "local_with_time_gain")
        self.assertGreater(row["FIN_equivalence_band_rate"], 0.95)

    def test_not_strict_validation(self):
        self.assertFalse(P["159"]["strict_FIN_validation"])


class TestProgram160(unittest.TestCase):
    def test_legacy_period(self):
        self.assertEqual(P["160"]["legacy_period"], 8)

    def test_strict_aperiodic(self):
        self.assertFalse(P["160"]["strict_has_nonzero_integer_period"])

    def test_not_z12_character(self):
        self.assertFalse(P["160"]["strict_is_Z12_character"])
        self.assertGreater(P["160"]["Z12_character_defect"], 1)

    def test_nonlinear_map_not_excluded(self):
        self.assertFalse(P["160"]["non_equivariant_or_affine_completion_excluded"])


class TestProgram161(unittest.TestCase):
    def test_grammar_exhausted(self):
        self.assertEqual(P["161"]["formulas_enumerated"], 343)

    def test_solutions_exist(self):
        self.assertGreater(P["161"]["number_of_solutions"], 0)

    def test_unique_shortest(self):
        self.assertTrue(P["161"]["unique_minimum_is_dimension_minus_one"])
        self.assertEqual(P["161"]["minimum_L1_complexity"], 1)

    def test_not_strict_source(self):
        self.assertFalse(P["161"]["strict_source_theorem"])


class TestProgram162(unittest.TestCase):
    def test_six_edges(self):
        self.assertEqual(len(P["162"]["edges"]), 6)

    def test_no_roles_transferable(self):
        self.assertEqual(P["162"]["eligible_physical_roles_count"], 0)
        self.assertTrue(not any(P["162"]["role_transfer_eligible"].values()))

    def test_matrix_has_expected_roles(self):
        self.assertIn("beta_tors^N hierarchy", P["162"]["roles"])


class TestProgram163(unittest.TestCase):
    def test_protocol_hash_matches_previous(self):
        previous = json.loads(
            (ROOT / "FIN_Programs_138_150_PreData_Physical_Protocol.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(
            P["163"]["protocol_core_sha256"], previous["canonical_core_sha256"]
        )

    def test_no_external_data(self):
        self.assertEqual(P["163"]["admitted_datasets"], 0)
        self.assertFalse(P["163"]["external_data_admitted"])

    def test_validator_accepts_complete_json(self):
        required = {"record_id", "time"}
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "candidate.json"
            path.write_text('[{"record_id": 1, "time": 1.0}]', encoding="utf-8")
            self.assertTrue(
                research.validate_external_candidate(path, required)["admitted"]
            )

    def test_validator_rejects_missing(self):
        required = {"record_id", "time"}
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "candidate.json"
            path.write_text('[{"record_id": 1}]', encoding="utf-8")
            result = research.validate_external_candidate(path, required)
            self.assertFalse(result["admitted"])
            self.assertEqual(result["missing"], ["time"])


class TestProgram164(unittest.TestCase):
    def test_six_axioms(self):
        self.assertEqual(set(P["164"]["axioms"]), {f"A{i}" for i in range(6)})

    def test_all_subsets(self):
        self.assertEqual(len(P["164"]["all_64_subsets"]), 64)

    def test_minimal_math_core(self):
        self.assertEqual(
            P["164"]["layered_minima"]["minimal_mathematical_core"], ["A0"]
        )

    def test_signed_operational_minimum(self):
        self.assertEqual(
            P["164"]["layered_minima"][
                "signed_dimensionless_operational_physics"
            ],
            ["A0", "A1", "A2", "A3"],
        )

    def test_calibrated_signed_minimum(self):
        self.assertEqual(
            P["164"]["layered_minima"]["calibrated_signed_test"],
            ["A0", "A1", "A2", "A3", "A4"],
        )

    def test_generic_operational_does_not_need_signed_resource(self):
        self.assertEqual(
            P["164"]["layered_minima"][
                "generic_dimensionless_operational_physics"
            ],
            ["A0", "A1", "A3"],
        )

    def test_every_axiom_has_countermodel(self):
        self.assertEqual(
            set(P["164"]["removal_countermodels"]), {f"A{i}" for i in range(6)}
        )

    def test_no_single_universal_theorem(self):
        self.assertFalse(
            P["164"]["smallest_single_theorem_connecting_every_layer_found"]
        )


class TestGlobalGuardrails(unittest.TestCase):
    def test_no_forbidden_promotions(self):
        verdict = RESULTS["global_verdict"]
        for key in [
            "strict_selector_closed",
            "QW_2191_discharged",
            "internal_units_derived",
            "legacy_strict_bridge_completed",
            "role_transfer_started",
            "L_total_derived",
            "ToE_claimed",
        ]:
            self.assertFalse(verdict[key])

    def test_author_spelling(self):
        self.assertEqual(RESULTS["metadata"]["creator"], "Żuchowski, Krzysztof")

    def test_no_firecrawl_or_external_data(self):
        self.assertFalse(RESULTS["metadata"]["firecrawl_used"])
        self.assertFalse(RESULTS["metadata"]["external_data_used"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
