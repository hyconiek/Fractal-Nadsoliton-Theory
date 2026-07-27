"""Regression and scientific-boundary tests for FIN Programs 178--190."""

from pathlib import Path
import json
import math
import unittest

import numpy as np

import fin_programs_178_190_composition_process_scale as study


ROOT = Path(__file__).resolve().parent
DATA = json.loads(
    (ROOT / "FIN_Programs_178_190_Composition_Process_Scale_Results.json").read_text(
        encoding="utf-8"
    )
)
P = DATA["programs"]


class TestGlobalGuardrails(unittest.TestCase):
    def test_author(self):
        self.assertEqual(DATA["metadata"]["creator"], "Żuchowski, Krzysztof")

    def test_no_firecrawl(self):
        self.assertFalse(DATA["metadata"]["firecrawl_used"])

    def test_no_external_web(self):
        self.assertFalse(DATA["metadata"]["external_web_used"])

    def test_no_selector_promotion(self):
        self.assertFalse(DATA["global_verdict"]["QW_2191_discharged"])

    def test_no_bridge_promotion(self):
        self.assertFalse(DATA["global_verdict"]["legacy_strict_bridge_completed"])

    def test_no_role_transfer(self):
        self.assertFalse(DATA["global_verdict"]["legacy_role_transfer_started"])

    def test_no_toe(self):
        self.assertFalse(DATA["global_verdict"]["L_total_or_ToE_claimed"])


class TestProgram178(unittest.TestCase):
    def setUp(self):
        self.d = P["178"]
        self.rows = {r["candidate"]: r for r in self.d["candidate_rows"]}

    def test_raw_cardinality_fails(self):
        self.assertFalse(
            self.rows["raw_cardinality_alpha_over_h"]["replication_intensive"]
        )

    def test_ratio_is_replication_intensive(self):
        self.assertTrue(
            self.rows["replication_invariant_alpha_over_log_h"][
                "replication_intensive"
            ]
        )

    def test_ratio_fails_label_split(self):
        self.assertFalse(
            self.rows["replication_invariant_alpha_over_log_h"][
                "label_split_invariant"
            ]
        )

    def test_constant_passes_both(self):
        row = self.rows["constant_log2"]
        self.assertTrue(row["replication_intensive"])
        self.assertTrue(row["label_split_invariant"])

    def test_no_unique_value(self):
        self.assertFalse(self.d["unique_target_free_beta_selected"])


class TestProgram179(unittest.TestCase):
    def setUp(self):
        self.d = P["179"]

    def test_all_positive_betas_collapse(self):
        self.assertLess(
            max(x["collapse_residual"] for x in self.d["orbit_collapse"]), 1e-12
        )

    def test_no_strict_source(self):
        self.assertEqual(self.d["accepted_strict_sources"], 0)

    def test_torsor_carrier(self):
        self.assertEqual(self.d["constructed_object"]["carrier"], "R_{>0}")

    def test_coordinate_choice_not_source(self):
        row = next(
            x for x in self.d["source_candidates"] if x["source"] == "coordinate_choice_beta=1"
        )
        self.assertIn("gauge", row["verdict"])


class TestProgram180(unittest.TestCase):
    def setUp(self):
        self.d = P["180"]

    def test_target_is_three_percent(self):
        self.assertAlmostEqual(self.d["preregistered_target"], 0.03)

    def test_certificate_open(self):
        self.assertFalse(self.d["full_one_engine_certificate"])

    def test_arb_absent(self):
        self.assertFalse(self.d["arb_or_python_flint_available"])

    def test_cancellation_small(self):
        self.assertLess(self.d["analytic_cancellation_bound_from_P165"], 1e-6)

    def test_formal_width_above_target(self):
        self.assertGreater(
            self.d["formal_worst_relative_enclosure_from_P151"],
            self.d["preregistered_target"],
        )


class TestProgram181(unittest.TestCase):
    def setUp(self):
        self.d = P["181"]

    def test_exact_identity(self):
        self.assertTrue(self.d["exact_identity"])

    def test_exact_energy(self):
        self.assertEqual(self.d["exact_quadratic_form"], "15")
        self.assertEqual(self.d["exact_dirichlet_sum"], "15")

    def test_constant_kernel(self):
        self.assertTrue(self.d["exact_constant_kernel"])

    def test_eigenvalues(self):
        self.assertTrue(
            np.allclose(self.d["eigenvalues"], [0, 1, 1, 2], atol=1e-12)
        )

    def test_unitarity(self):
        self.assertLess(self.d["max_unitarity_residual"], 2e-14)

    def test_stochasticity(self):
        self.assertLess(self.d["max_stochastic_row_residual"], 2e-14)
        self.assertGreaterEqual(self.d["minimum_heat_entry"], -2e-15)

    def test_no_false_lean_claim(self):
        self.assertFalse(self.d["lean_machine_compiled"])
        self.assertTrue((ROOT / self.d["lean_source"]).exists())


class TestProgram182(unittest.TestCase):
    def setUp(self):
        self.d = P["182"]

    def test_counterexample_same_M(self):
        self.assertTrue(self.d["counterexample"]["same_M"])

    def test_counterexample_infeasible(self):
        self.assertFalse(self.d["counterexample"]["conversion_feasible"])

    def test_reachable_M_bound(self):
        self.assertAlmostEqual(
            self.d["counterexample"]["maximum_reachable_target_M"], 0.36, places=6
        )

    def test_M_false_positives_exist(self):
        self.assertGreater(self.d["grid_false_positives_of_M_only"], 0)

    def test_direct_formula(self):
        factor, _, _ = study.max_covariant_transverse_factor(0.0, 0.8)
        self.assertAlmostEqual(factor, 0.6, places=5)


class TestProgram183(unittest.TestCase):
    def setUp(self):
        self.d = P["183"]

    def test_effective_sample_size(self):
        self.assertEqual(
            self.d["nominal_n"],
            self.d["independent_blocks"] * self.d["block_size"],
        )

    def test_valid_epsilon_larger(self):
        self.assertGreater(
            self.d["epsilon_valid_block"], self.d["epsilon_invalid_nominal"]
        )

    def test_valid_coverage(self):
        self.assertGreaterEqual(self.d["coverage_valid_block"], 0.95)

    def test_nominal_undercoverage(self):
        self.assertLess(self.d["coverage_invalid_nominal"], 0.85)

    def test_run_diagnostic(self):
        self.assertEqual(
            self.d["run_length_detected_block_size"], self.d["block_size"]
        )


class TestProgram184(unittest.TestCase):
    def setUp(self):
        self.d = P["184"]
        self.s = self.d["simulation"]

    def test_full_rank(self):
        self.assertEqual(self.d["design_rank"], self.d["parameters"])

    def test_one_control_is_algebraically_sufficient_in_model(self):
        self.assertEqual(self.d["one_control_reduced_rank"], 5)

    def test_naive_confounded(self):
        self.assertGreater(abs(self.s["naive_mean"] - self.s["true_target"]), 0.5)

    def test_controlled_unbiased(self):
        self.assertLess(
            abs(self.s["controlled_mean"] - self.s["true_target"]), 0.01
        )

    def test_violation_diagnostic_power(self):
        self.assertGreater(
            self.s["shared_gain_violation_detection_rate"], 0.95
        )


class TestProgram185(unittest.TestCase):
    def setUp(self):
        self.d = P["185"]
        self.rows = {r["unknown_model"]: r for r in self.d["rows"]}

    def test_protocol_exists(self):
        self.assertTrue((ROOT / self.d["protocol"]).exists())

    def test_exact_order_blindness(self):
        self.assertEqual(self.d["exact_order_blindness_residual"], 0.0)

    def test_ordered_unknown_false_accepts(self):
        self.assertGreater(
            self.rows["copula_ordered_gaussian"]["false_accept_rate"], 0.85
        )

    def test_open_set_not_complete(self):
        self.assertFalse(self.d["open_set_complete"])

    def test_no_external_data(self):
        self.assertFalse(self.d["external_data_used"])


class TestProgram186(unittest.TestCase):
    def setUp(self):
        self.d = P["186"]

    def test_memory_witness_nonzero(self):
        self.assertGreater(abs(self.d["nonfactorization_memory_witness"]), 0.1)

    def test_static_echo(self):
        self.assertAlmostEqual(self.d["echo_static"], 1.0)

    def test_markov_echo_not_restored(self):
        self.assertLess(self.d["echo_markov"], 1.0)

    def test_memory_flag(self):
        self.assertTrue(self.d["memory_present_for_quasistatic_process"])

    def test_four_intervention_vectors(self):
        self.assertEqual(len(self.d["intervention_vectors"]), 4)


class TestProgram187(unittest.TestCase):
    def setUp(self):
        self.d = P["187"]

    def test_same_A_distinct_channels(self):
        self.assertTrue(self.d["same_A_distinct_channels"])

    def test_endpoints(self):
        rows = self.d["sample_rows"]
        self.assertAlmostEqual(rows[0]["gamma"], 1.0)
        self.assertAlmostEqual(rows[-1]["gamma"], 0.0, places=12)

    def test_choi_positive(self):
        for row in self.d["sample_rows"]:
            self.assertGreaterEqual(min(row["choi_eigenvalues"]), -1e-12)

    def test_operational_replacement(self):
        self.assertIn("process tensor", self.d["minimal_operational_replacement"])


class TestProgram188(unittest.TestCase):
    def setUp(self):
        self.d = P["188"]

    def test_residual(self):
        self.assertLess(self.d["correction_residual"], 1e-14)

    def test_no_phase_source(self):
        self.assertFalse(self.d["strict_phase_source_exported"])

    def test_delta_frequency(self):
        self.assertAlmostEqual(
            self.d["delta_omega"], 0.18575 - math.pi / 4
        )

    def test_cohomology(self):
        self.assertIn("H^1", self.d["cohomology"])


class TestProgram189(unittest.TestCase):
    def setUp(self):
        self.d = P["189"]

    def test_unitary_invariant(self):
        self.assertLess(self.d["max_unitary_orbit_residual"], 2e-14)

    def test_heat_invariant(self):
        self.assertLess(self.d["max_heat_orbit_residual"], 2e-14)

    def test_resolvent_covariant(self):
        self.assertLess(self.d["max_resolvent_covariance_residual"], 2e-12)

    def test_normalized_gap_invariant(self):
        self.assertLess(abs(self.d["normalized_gap_variation"]), 1e-12)

    def test_raw_gap_changes(self):
        self.assertGreater(self.d["raw_gap_scale_ratio"], 1e5)

    def test_clock_is_conditional(self):
        self.assertIn("conditional", self.d["constructed_object"]["status"])


class TestProgram190(unittest.TestCase):
    def setUp(self):
        self.d = P["190"]

    def test_three_bundles_audited(self):
        self.assertEqual(len(self.d["audited_bundles"]), 3)

    def test_no_admission(self):
        self.assertEqual(self.d["admitted_bundles"], [])

    def test_intake_failure_record(self):
        self.assertTrue(self.d["intake_failure_record"])

    def test_no_external_execution(self):
        self.assertFalse(self.d["external_protocol_executed"])

    def test_reference_control_required(self):
        self.assertIn("reference_control", self.d["required_fields"])


if __name__ == "__main__":
    unittest.main()
