"""Regression and scientific-boundary tests for FIN Programs 191--203."""

from fractions import Fraction
from pathlib import Path
import csv
import hashlib
import json
import math
import unittest

import numpy as np

import fin_programs_191_203_reference_process_prediction as study


ROOT = Path(__file__).resolve().parent
DATA = json.loads(
    (ROOT / "FIN_Programs_191_203_Reference_Process_Prediction_Results.json").read_text(
        encoding="utf-8"
    )
)
P = DATA["programs"]


class TestGlobalGuardrails(unittest.TestCase):
    def test_all_programs_executed(self):
        self.assertEqual(
            DATA["global_verdict"]["programs_executed"], list(range(191, 204))
        )

    def test_creator_spelling(self):
        self.assertEqual(DATA["metadata"]["creator"], "Żuchowski, Krzysztof")

    def test_no_firecrawl(self):
        self.assertFalse(DATA["metadata"]["firecrawl_used"])

    def test_no_external_web(self):
        self.assertFalse(DATA["metadata"]["external_web_used"])

    def test_no_external_data(self):
        self.assertFalse(DATA["metadata"]["external_dataset_used"])

    def test_no_selector_promotion(self):
        self.assertFalse(DATA["global_verdict"]["QW_2191_discharged"])
        self.assertFalse(DATA["global_verdict"]["strict_selector_exported"])

    def test_no_unit_promotion(self):
        self.assertFalse(DATA["global_verdict"]["canonical_physical_unit_exported"])

    def test_no_bridge_or_role_transfer(self):
        self.assertFalse(DATA["global_verdict"]["legacy_strict_bridge_completed"])
        self.assertFalse(DATA["global_verdict"]["legacy_role_transfer_started"])

    def test_no_toe_or_external_validation(self):
        self.assertFalse(DATA["global_verdict"]["L_total_or_ToE_claimed"])
        self.assertFalse(
            DATA["global_verdict"]["external_physical_validation_claimed"]
        )


class TestProgram191(unittest.TestCase):
    def setUp(self):
        self.d = P["191"]

    def test_exact_expectations(self):
        self.assertAlmostEqual(self.d["uniform_central_expectation"], 9 / 5)
        self.assertAlmostEqual(self.d["normalized_full_trace_expectation"], 17 / 9)

    def test_weights_normalized(self):
        self.assertAlmostEqual(sum(self.d["uniform_central_weights"]), 1)
        self.assertAlmostEqual(sum(self.d["normalized_full_trace_weights"]), 1)

    def test_isomorphic_blocks_invariant(self):
        self.assertEqual(self.d["permutation_residual"], 0)

    def test_tensor_trace(self):
        self.assertLess(self.d["simple_trace_tensor_residual"], 1e-15)

    def test_no_unique_state(self):
        self.assertFalse(self.d["unique_reference_state_selected"])


class TestProgram192(unittest.TestCase):
    def setUp(self):
        self.d = P["192"]

    def test_orbit_collapse(self):
        self.assertLess(self.d["maximum_orbit_residual"], 1e-12)

    def test_invariant_frequency(self):
        strict = self.d["invariant_rows"]["strict"]
        legacy = self.d["invariant_rows"]["legacy"]
        self.assertAlmostEqual(strict["nu"], study.OMEGA_STRICT)
        self.assertAlmostEqual(
            legacy["nu"], study.OMEGA_LEGACY / study.BETA_LEGACY
        )

    def test_shape_gap_nonzero(self):
        self.assertGreater(self.d["quotient_profile_L2_gap"], 0.3)

    def test_beta_one_only_gauge(self):
        self.assertTrue(self.d["beta_one_is_gauge_representative"])
        self.assertFalse(self.d["gauge_closes_legacy_strict_bridge"])


class TestProgram193(unittest.TestCase):
    def setUp(self):
        self.d = P["193"]

    def test_dependency_count(self):
        self.assertEqual(self.d["total_components"], 5)
        self.assertLess(self.d["closed_components"], self.d["total_components"])

    def test_no_false_certificate(self):
        self.assertFalse(self.d["full_one_engine_certificate"])
        self.assertFalse(self.d["formal_sub_three_percent_certificate"])

    def test_existing_width_above_target(self):
        self.assertGreater(
            self.d["formal_worst_relative_enclosure"],
            self.d["target_relative_width"],
        )

    def test_environment_explicit(self):
        self.assertIn("python_flint_available", self.d["environment"])


class TestProgram194(unittest.TestCase):
    def setUp(self):
        self.d = P["194"]

    def test_many_exact_cases(self):
        self.assertEqual(self.d["exact_rational_cycle_cases"], 150)

    def test_fraction_identity_direct(self):
        w = study._cycle_matrix_fraction(5)
        f = [Fraction(i - 2, 3) for i in range(5)]
        q, rhs, _ = study._dirichlet_exact(w, f)
        self.assertEqual(q, rhs)
        self.assertGreaterEqual(q, 0)

    def test_constant_kernel_direct(self):
        w = study._cycle_matrix_fraction(7)
        f = [Fraction(5, 3)] * 7
        q, rhs, af = study._dirichlet_exact(w, f)
        self.assertEqual(q, 0)
        self.assertEqual(rhs, 0)
        self.assertTrue(all(x == 0 for x in af))

    def test_lean_status_honest(self):
        self.assertFalse(self.d["lean_machine_compiled"])
        self.assertTrue((ROOT / self.d["lean_source"]).exists())

    def test_numerical_semigroups(self):
        self.assertLess(self.d["maximum_unitarity_residual"], 1e-12)
        self.assertLess(self.d["maximum_heat_row_residual"], 1e-12)
        self.assertGreater(self.d["minimum_heat_entry"], -1e-12)


class TestProgram195(unittest.TestCase):
    def setUp(self):
        self.d = P["195"]

    def test_counterexample_boundary(self):
        self.assertAlmostEqual(self.d["counterexample_boundary"], 0.36, places=12)
        self.assertEqual(self.d["counterexample_active_constraint"], "source-kink")

    def test_closed_formula_matches_old_grid(self):
        self.assertLess(self.d["maximum_closed_vs_choi_grid_residual"], 1e-5)

    def test_au_margin_on_boundary(self):
        boundary, _, _ = study.reflection_boundary_closed_form(0.6, 0.0, 0.8)
        margin = study.alberti_uhlmann_margin((0.6, 0.0), (boundary, 0.8))
        self.assertGreaterEqual(margin, -1e-8)

    def test_above_boundary_fails(self):
        margin = study.alberti_uhlmann_margin((0.6, 0.0), (0.37, 0.8))
        self.assertLess(margin, -1e-4)

    def test_tensor_catalysis_not_claimed(self):
        self.assertFalse(self.d["tensor_catalysis_classified"])


class TestProgram196(unittest.TestCase):
    def setUp(self):
        self.d = P["196"]

    def test_mixing_bound_contracts(self):
        p = self.d["refresh_probability"]
        self.assertLess((1 - p) ** 100, 0.01)

    def test_retained_fraction(self):
        self.assertAlmostEqual(
            self.d["mean_retained_fraction"],
            self.d["refresh_probability"],
            delta=0.002,
        )

    def test_valid_coverage(self):
        self.assertGreaterEqual(self.d["regenerative_coverage"], 0.95)

    def test_nominal_undercoverage(self):
        self.assertLess(self.d["nominal_record_coverage"], 0.95)

    def test_valid_interval_cost(self):
        self.assertGreater(
            self.d["regenerative_mean_width"], self.d["nominal_mean_width"]
        )


class TestProgram197(unittest.TestCase):
    def setUp(self):
        self.d = P["197"]
        self.rows = {row["challenge"]: row for row in self.d["rows"]}

    def test_protocol_hash(self):
        path = ROOT / self.d["protocol"]
        self.assertTrue(path.exists())
        record = json.loads(path.read_text(encoding="utf-8"))
        self.assertEqual(
            study.canonical_digest(record["core"]), record["canonical_core_sha256"]
        )

    def test_multiset_identity(self):
        self.assertEqual(self.d["exact_multiset_residual_sorted_challenge"], 0)

    def test_iid_specificity(self):
        self.assertLess(self.rows["iid_gaussian_validation"]["rejection_rate"], 0.05)
        self.assertLess(self.rows["iid_stable_0.8"]["rejection_rate"], 0.05)

    def test_temporal_alternatives_rejected(self):
        for name in [
            "sorted_same_multiset",
            "AR1_phi_0.8",
            "repeated_block_20",
            "linear_drift",
        ]:
            self.assertGreater(self.rows[name]["rejection_rate"], 0.95)

    def test_no_completeness(self):
        self.assertFalse(self.d["closed_set_completeness_claimed"])


class TestProgram198(unittest.TestCase):
    def setUp(self):
        self.d = P["198"]

    def test_minimal_rank(self):
        self.assertEqual(self.d["rank_three_contrasts"], 3)
        self.assertEqual(self.d["rank_without_echo"], 2)

    def test_determinant(self):
        matrix = np.asarray(self.d["constructed_object"]["linear_log_model"])
        self.assertAlmostEqual(abs(np.linalg.det(matrix)), 1)

    def test_estimates_unbiased(self):
        sim = self.d["simulation"]
        self.assertTrue(
            np.allclose(sim["mean_estimate"], sim["truth"], atol=0.002)
        )

    def test_condition_finite(self):
        self.assertLess(self.d["condition_number"], 10)


class TestProgram199(unittest.TestCase):
    def setUp(self):
        self.d = P["199"]
        self.c = self.d["one_time_counterexample"]

    def test_one_time_same(self):
        self.assertAlmostEqual(
            self.c.get("one_time_characteristic_A", 0.8),
            self.c.get("one_time_characteristic_B", 0.8),
        )

    def test_two_time_different(self):
        self.assertGreater(
            abs(
                self.c["two_time_characteristic_A"]
                - self.c["two_time_characteristic_B"]
            ),
            0.7,
        )

    def test_process_distinguishable(self):
        self.assertTrue(self.d["multi_time_processes_distinguishable"])


class TestProgram200(unittest.TestCase):
    def setUp(self):
        self.d = P["200"]

    def test_no_accepted_source(self):
        self.assertEqual(self.d["accepted_source_operations"], 0)
        self.assertFalse(self.d["strict_phase_source_exported"])

    def test_target_coded_candidates_flagged(self):
        coded = [
            row for row in self.d["candidate_rows"] if row["produces_strict"]
        ]
        self.assertTrue(coded)
        self.assertTrue(all(not row["target_independent"] for row in coded))

    def test_endomorphism_images_miss_target(self):
        self.assertGreater(self.d["minimum_distance_from_endomorphism_images"], 0.1)


class TestProgram201(unittest.TestCase):
    def setUp(self):
        self.d = P["201"]

    def test_invariants_constant(self):
        self.assertLess(max(self.d["maximum_invariance_residuals"].values()), 1e-12)

    def test_raw_gap_changes(self):
        self.assertGreater(self.d["raw_gap_ratio"], 1e11)

    def test_physical_scale_not_exported(self):
        self.assertFalse(self.d["strict_physical_scale_exported"])

    def test_raw_quantities_excluded(self):
        self.assertIn("raw spectral gap", self.d["noninvariants"])


class TestProgram202(unittest.TestCase):
    def setUp(self):
        self.d = P["202"]

    def test_ten_of_eleven(self):
        self.assertEqual(self.d["passed_required_fields"], 10)
        self.assertEqual(self.d["total_required_fields"], 11)
        self.assertEqual(
            self.d["failed_fields"], ["independent_analysis_boundary"]
        )

    def test_external_not_admitted(self):
        self.assertFalse(self.d["external_bundle_admitted"])
        self.assertFalse(self.d["external_protocol_executed"])

    def test_bundle_hash(self):
        path = ROOT / self.d["constructed_object"]["raw_event_file"]
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        self.assertEqual(
            digest, self.d["constructed_object"]["raw_event_sha256"]
        )

    def test_event_order(self):
        path = ROOT / self.d["constructed_object"]["raw_event_file"]
        with path.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        sequence = [int(row["sequence"]) for row in rows]
        self.assertEqual(sequence, list(range(len(sequence))))
        self.assertGreater(len(rows), 10000)


class TestProgram203(unittest.TestCase):
    def setUp(self):
        self.d = P["203"]

    def test_gate_remains_failed(self):
        self.assertFalse(self.d["program202_external_gate_passed"])
        self.assertFalse(self.d["external_physical_prediction_tested"])

    def test_dry_run_only(self):
        self.assertTrue(self.d["synthetic_method_validation_only"])

    def test_heldout_prediction_reasonable(self):
        self.assertLess(abs(self.d["heldout_residual"]), 0.06)

    def test_estimated_parameters_physical_in_declared_model(self):
        p = self.d["estimated_parameters"]
        self.assertGreater(p["blur"], 0)
        self.assertLessEqual(p["blur"], 1)
        self.assertGreater(p["variance"], 0)
        self.assertLess(abs(p["covariance"]), p["variance"])

    def test_dependency_ledger_complete(self):
        self.assertEqual(set(self.d["dependency_ledger"]), {"W0", "CA", "OP"})


if __name__ == "__main__":
    unittest.main()
