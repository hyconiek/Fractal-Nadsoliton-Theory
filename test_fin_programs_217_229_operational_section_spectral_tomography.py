#!/usr/bin/env python3
"""Regression and scientific-boundary tests for FIN Programs 217--229."""

from fractions import Fraction
from pathlib import Path
import hashlib
import json
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "FIN_Programs_217_229_Operational_Section_Spectral_Tomography_Results.json"
DATA = json.loads(RESULTS.read_text(encoding="utf-8"))
P = DATA["programs"]
G = DATA["global_verdict"]
FIG = ROOT / "FIN_Programs_217_229_Operational_Section_Spectral_Tomography_Figures"


def verify_file_hash(record: dict) -> None:
    path = ROOT / record["path"]
    assert path.is_file()
    assert hashlib.sha256(path.read_bytes()).hexdigest() == record["sha256"]


def verify_canonical(path: Path) -> None:
    record = json.loads(path.read_text(encoding="utf-8"))
    payload = json.dumps(
        record["core"],
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
    ).encode("utf-8")
    assert hashlib.sha256(payload).hexdigest() == record["canonical_core_sha256"]


class TestGlobalGuardrails(unittest.TestCase):
    def test_program_range(self):
        self.assertEqual(G["programs_executed"], list(range(217, 230)))

    def test_thirteen_objects(self):
        self.assertEqual(len(G["new_theoretical_objects"]), 13)

    def test_creator(self):
        self.assertEqual(DATA["metadata"]["creator"], "Żuchowski, Krzysztof")

    def test_no_firecrawl_or_web(self):
        self.assertFalse(DATA["metadata"]["firecrawl_used"])
        self.assertFalse(DATA["metadata"]["external_web_used"])

    def test_no_external_dataset(self):
        self.assertFalse(DATA["metadata"]["external_dataset_used"])

    def test_no_selector(self):
        self.assertFalse(G["QW_2191_discharged"])
        self.assertFalse(G["strict_selector_exported"])

    def test_no_unit(self):
        self.assertFalse(G["canonical_physical_unit_exported"])

    def test_no_bridge_or_role_transfer(self):
        self.assertFalse(G["legacy_strict_bridge_completed"])
        self.assertFalse(G["legacy_role_transfer_started"])

    def test_no_toe(self):
        self.assertFalse(G["L_total_or_ToE_claimed"])

    def test_no_external_validation(self):
        self.assertFalse(G["external_physical_validation_claimed"])

    def test_all_figures(self):
        names = sorted(FIG.glob("program*.png"))
        self.assertEqual(len(names), 13)
        self.assertTrue(all(path.stat().st_size > 20_000 for path in names))


class TestProgram217(unittest.TestCase):
    def setUp(self):
        self.d = P["217"]

    def test_bistochastic_invariance(self):
        self.assertLess(self.d["maximum_uniform_invariance_residual"], 1e-14)

    def test_exact_uniform_expectation(self):
        self.assertEqual(self.d["uniform_dimension_expectation_exact"], "9/5")
        self.assertAlmostEqual(self.d["uniform_dimension_expectation"], 1.8)

    def test_prepared_state_differs(self):
        self.assertAlmostEqual(self.d["prepared_dimension_expectation"], 1.9)

    def test_hoeffding_coverage(self):
        self.assertGreaterEqual(self.d["empirical_simultaneous_coverage"], 0.95)

    def test_not_strict_source(self):
        self.assertFalse(self.d["strict_central_state_source_exported"])


class TestProgram218(unittest.TestCase):
    def setUp(self):
        self.d = P["218"]

    def test_polynomial_rank(self):
        self.assertEqual(self.d["polynomial_constraint_rank"], 8)
        self.assertEqual(self.d["polynomial_constraint_nullity"], 1)

    def test_only_linear_survives(self):
        self.assertEqual(self.d["surviving_monomial_degrees"], [1])

    def test_increment(self):
        self.assertEqual(self.d["strict_increment_exact"], "4/5")

    def test_only_zero_slope_fixes_target(self):
        self.assertEqual(self.d["sampled_slopes_fixing_target"], 1)

    def test_no_eta_source(self):
        self.assertFalse(self.d["target_independent_eta_source_exported"])


class TestProgram219(unittest.TestCase):
    def setUp(self):
        self.d = P["219"]

    def test_attestation_hash(self):
        verify_file_hash(self.d["constructed_object"])

    def test_attestation_canonical(self):
        verify_canonical(ROOT / self.d["constructed_object"]["path"])

    def test_engine_absent(self):
        self.assertFalse(self.d["python_flint_available"])
        self.assertFalse(self.d["arb_binary_available"])
        self.assertFalse(self.d["docker_server_accessible"])

    def test_no_false_certificate(self):
        self.assertFalse(self.d["formal_run_executed"])
        self.assertFalse(self.d["formal_sub_three_percent_certificate"])


class TestProgram220(unittest.TestCase):
    def setUp(self):
        self.d = P["220"]

    def test_record_hash_and_canonical(self):
        verify_file_hash(self.d["constructed_object"])
        verify_canonical(ROOT / self.d["constructed_object"]["path"])

    def test_source_hash(self):
        path = ROOT / self.d["dependency_free_source"]
        self.assertEqual(
            hashlib.sha256(path.read_bytes()).hexdigest(),
            self.d["dependency_free_source_sha256"],
        )

    def test_core_compiled(self):
        self.assertTrue(self.d["dependency_free_machine_compiled"])
        self.assertIn("30", self.d["dependency_free_output_tail"])

    def test_five_scoped_statements(self):
        self.assertEqual(len(self.d["machine_checked_statements"]), 5)

    def test_general_not_promoted(self):
        self.assertFalse(self.d["general_mathlib_library_compiled"])
        self.assertIn("Mathlib", self.d["general_output_tail"])


class TestProgram221(unittest.TestCase):
    def setUp(self):
        self.d = P["221"]

    def test_violation_all_imperfect_grid_points(self):
        self.assertEqual(self.d["positive_uncorrelated_violation_points"], 400)

    def test_perfect_boundary(self):
        self.assertAlmostEqual(self.d["perfect_reference_violation"], 0.0)

    def test_correlated_marginals(self):
        example = self.d["correlated_example"]
        self.assertAlmostEqual(example["returned_reference_marginal"], 0.8)
        self.assertAlmostEqual(example["target_marginal"], 0.6)

    def test_correlation_cost_positive(self):
        self.assertGreater(self.d["correlated_example"]["mutual_information_bits"], 0)

    def test_conditional_radius_physical(self):
        self.assertLessEqual(
            self.d["correlated_example"]["conditional_preparation_radius"], 1
        )

    def test_no_selector_source(self):
        self.assertFalse(self.d["strict_selector_source_exported"])


class TestProgram222(unittest.TestCase):
    def setUp(self):
        self.d = P["222"]

    def test_lower_bound_below_truth(self):
        self.assertLess(
            self.d["one_sided_calibrated_p_lower"],
            self.d["true_refresh_rate_for_challenge"],
        )

    def test_calibration_coverage(self):
        self.assertGreaterEqual(self.d["empirical_lower_bound_coverage"], 0.98)

    def test_budget_accounting(self):
        row = self.d["robust_design"]
        self.assertAlmostEqual(
            row["coupling_budget"] + row["dkw_budget"],
            self.d["per_sequence_failure_budget"],
        )

    def test_robust_design_improves_inherited(self):
        self.assertLess(
            self.d["robust_design"]["epsilon"],
            self.d["inherited_equal_split_epsilon"],
        )

    def test_oracle_is_better(self):
        self.assertLess(
            self.d["oracle_design"]["epsilon"],
            self.d["robust_design"]["epsilon"],
        )

    def test_nontrivial_thinning(self):
        self.assertGreater(self.d["robust_design"]["lag"], 1)
        self.assertLess(self.d["robust_design"]["count"], self.d["record_length"])


class TestProgram223(unittest.TestCase):
    def setUp(self):
        self.d = P["223"]

    def test_time_uniform_null(self):
        self.assertLessEqual(self.d["online_rank_null_rejection_rate"], 0.05)

    def test_alternative_power(self):
        self.assertGreaterEqual(self.d["online_rank_shifted_rejection_rate"], 0.95)

    def test_fast_alternative_crossing(self):
        self.assertLess(self.d["online_rank_shifted_mean_crossing_time"], 5)

    def test_naive_conditional_failure_has_positive_measure(self):
        self.assertGreater(
            self.d["naive_fraction_with_conditional_mean_e_above_one"], 0.25
        )

    def test_naive_conditional_failure_can_be_large(self):
        self.assertGreater(self.d["naive_maximum_conditional_mean_e"], 1.5)

    def test_rank_filtration_declared(self):
        self.assertIn("filtration", self.d["theorem"])


class TestProgram224(unittest.TestCase):
    def setUp(self):
        self.d = P["224"]

    def test_coverage_near_nominal(self):
        self.assertGreater(self.d["simultaneous_empirical_coverage"], 0.93)
        self.assertLess(self.d["simultaneous_empirical_coverage"], 0.99)

    def test_width_reduction(self):
        for key, value in self.d["mean_parameter_interval_widths"].items():
            self.assertLess(value, self.d["inherited_exact_region_mean_widths"][key])

    def test_memory_power(self):
        self.assertGreaterEqual(self.d["memory_detection_power"], 0.95)

    def test_control_size(self):
        self.assertGreater(self.d["independent_control_null_rejection_rate"], 0.015)
        self.assertLess(self.d["independent_control_null_rejection_rate"], 0.085)

    def test_misspecification_power(self):
        self.assertGreater(
            self.d["independent_control_misspecified_rejection_rate"], 0.30
        )

    def test_not_exact_claim(self):
        self.assertIn("not an exact", self.d["claim_boundary"])


class TestProgram225(unittest.TestCase):
    def setUp(self):
        self.d = P["225"]

    def test_exact_interval(self):
        self.assertEqual(self.d["sharp_c3_lower_exact"], "-11/180")
        self.assertEqual(self.d["sharp_c3_upper_exact"], "11/20")

    def test_ordering(self):
        self.assertLess(self.d["sharp_c3_lower"], self.d["sharp_c3_upper"])

    def test_lower_witness(self):
        self.assertEqual(self.d["lower_witness_c3_exact"], "-11/180")

    def test_upper_witness(self):
        self.assertEqual(self.d["upper_witness_c3_exact"], "11/20")

    def test_lower_weights_sum(self):
        self.assertEqual(Fraction(11, 335) + 2 * Fraction(162, 335), 1)

    def test_upper_weights_sum(self):
        self.assertEqual(Fraction(11, 15) + 2 * Fraction(2, 15), 1)

    def test_witness_first_two_x_moments(self):
        self.assertEqual(self.d["lower_witness_moments_x"][:2], ["4/5", "3/4"])
        self.assertEqual(self.d["upper_witness_moments_x"][:2], ["4/5", "3/4"])


class TestProgram226(unittest.TestCase):
    def setUp(self):
        self.d = P["226"]

    def test_record_hash_and_canonical(self):
        verify_file_hash(self.d["constructed_object"])
        verify_canonical(ROOT / self.d["constructed_object"]["path"])

    def test_source_hash(self):
        path = ROOT / self.d["lean_source"]
        self.assertEqual(
            hashlib.sha256(path.read_bytes()).hexdigest(),
            self.d["lean_source_sha256"],
        )

    def test_finite_core_compiled(self):
        self.assertTrue(self.d["finite_torsion_core_machine_compiled"])
        self.assertIn("[0, 1, 2, 3, 4, 5, 6, 7]", self.d["finite_output_tail"])

    def test_analytic_scope_not_promoted(self):
        self.assertFalse(self.d["analytic_automatic_continuity_machine_formalized"])
        self.assertFalse(self.d["transcendence_non_torsion_machine_formalized"])

    def test_no_phase_source(self):
        self.assertFalse(self.d["strict_phase_source_exported"])


class TestProgram227(unittest.TestCase):
    def setUp(self):
        self.d = P["227"]

    def test_exact_reconstruction(self):
        self.assertLess(self.d["exact_reconstruction_signature_residual"], 1e-12)

    def test_positive_transition(self):
        self.assertGreater(self.d["transition_minimum_eigenvalue"], 0)

    def test_no_spd_failures(self):
        self.assertTrue(all(row["spd_failures"] == 0 for row in self.d["shot_rows"]))

    def test_accuracy_improves(self):
        means = [row["mean_fingerprint_distance"] for row in self.d["shot_rows"]]
        self.assertTrue(all(a > b for a, b in zip(means, means[1:])))

    def test_high_shot_pass(self):
        self.assertGreaterEqual(
            self.d["shot_rows"][-1]["program214_acceptance_rate"], 0.95
        )

    def test_low_shot_not_falsely_strong(self):
        self.assertLessEqual(
            self.d["shot_rows"][0]["program214_acceptance_rate"], 0.20
        )

    def test_no_external_data(self):
        self.assertFalse(self.d["external_data_used"])


class TestProgram228(unittest.TestCase):
    def setUp(self):
        self.d = P["228"]

    def test_request_hash_and_canonical(self):
        verify_file_hash(self.d["constructed_object"])
        verify_canonical(ROOT / self.d["constructed_object"]["path"])

    def test_no_admitted_bundle(self):
        self.assertEqual(self.d["admitted_bundles"], [])
        self.assertFalse(self.d["external_11_of_11_bundle_found"])

    def test_no_double_slit_record(self):
        self.assertFalse(self.d["double_slit_event_bundle_found"])

    def test_no_web(self):
        self.assertFalse(self.d["web_or_firecrawl_used"])

    def test_synthetic_not_admitted(self):
        synthetic = [row for row in self.d["candidate_rows"] if "Synthetic" in row["bundle"]][0]
        self.assertFalse(synthetic["admitted"])


class TestProgram229(unittest.TestCase):
    def setUp(self):
        self.d = P["229"]

    def test_lock_hash_and_canonical(self):
        verify_file_hash(self.d["constructed_object"])
        verify_canonical(ROOT / self.d["constructed_object"]["path"])

    def test_gate_failed(self):
        self.assertFalse(self.d["program228_gate_passed"])

    def test_not_executed(self):
        self.assertFalse(self.d["external_prediction_executed"])

    def test_no_validation(self):
        self.assertFalse(self.d["external_physical_validation_claimed"])


if __name__ == "__main__":
    unittest.main()
