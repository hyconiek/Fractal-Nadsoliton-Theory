"""Regression and claim-boundary tests for FIN Programs 204--216."""

from fractions import Fraction
from pathlib import Path
import hashlib
import json
import math
import unittest

import numpy as np

import fin_programs_204_216_categorical_catalytic_external as study


ROOT = Path(__file__).resolve().parent
DATA = json.loads(
    (ROOT / "FIN_Programs_204_216_Categorical_Catalytic_External_Results.json").read_text(
        encoding="utf-8"
    )
)
P = DATA["programs"]


def verify_canonical_record(path: Path) -> None:
    record = json.loads(path.read_text(encoding="utf-8"))
    assert study.canonical_digest(record["core"]) == record["canonical_core_sha256"]


class TestGlobalGuardrails(unittest.TestCase):
    def test_program_range(self):
        self.assertEqual(
            DATA["global_verdict"]["programs_executed"], list(range(204, 217))
        )

    def test_creator(self):
        self.assertEqual(DATA["metadata"]["creator"], "Żuchowski, Krzysztof")

    def test_no_web_or_firecrawl(self):
        self.assertFalse(DATA["metadata"]["firecrawl_used"])
        self.assertFalse(DATA["metadata"]["external_web_used"])

    def test_no_external_dataset(self):
        self.assertFalse(DATA["metadata"]["external_dataset_used"])

    def test_no_selector(self):
        self.assertFalse(DATA["global_verdict"]["QW_2191_discharged"])
        self.assertFalse(DATA["global_verdict"]["strict_selector_exported"])

    def test_no_unit(self):
        self.assertFalse(DATA["global_verdict"]["canonical_physical_unit_exported"])

    def test_no_bridge(self):
        self.assertFalse(DATA["global_verdict"]["legacy_strict_bridge_completed"])
        self.assertFalse(DATA["global_verdict"]["legacy_role_transfer_started"])

    def test_no_toe(self):
        self.assertFalse(DATA["global_verdict"]["L_total_or_ToE_claimed"])

    def test_no_external_validation(self):
        self.assertFalse(
            DATA["global_verdict"]["external_physical_validation_claimed"]
        )

    def test_thirteen_objects(self):
        self.assertEqual(len(DATA["global_verdict"]["new_theoretical_objects"]), 13)


class TestProgram204(unittest.TestCase):
    def setUp(self):
        self.d = P["204"]

    def test_rank_and_nullity(self):
        self.assertEqual(self.d["morita_constraint_rank"], 4)
        self.assertEqual(self.d["morita_constraint_nullity"], 1)

    def test_uniform_morita_weights(self):
        self.assertTrue(
            np.allclose(self.d["normalized_unique_weights"], np.full(5, 0.2))
        )

    def test_exact_expectations(self):
        self.assertAlmostEqual(self.d["morita_expectation"], 9 / 5)
        self.assertAlmostEqual(self.d["represented_trace_expectation"], 17 / 9)

    def test_all_hom_no_go(self):
        self.assertAlmostEqual(
            self.d["all_unital_homomorphism_no_go"]["defect"], 1 / 6
        )

    def test_not_strict_source(self):
        self.assertFalse(self.d["strict_central_measure_source_exported"])


class TestProgram205(unittest.TestCase):
    def setUp(self):
        self.d = P["205"]

    def test_increment(self):
        self.assertEqual(self.d["required_increment_exact"], "4/5")
        self.assertAlmostEqual(self.d["strict_eta"] - self.d["legacy_eta"], 0.8)

    def test_cocycle_residual(self):
        self.assertLess(self.d["maximum_cocycle_residual"], 1e-12)

    def test_no_candidate_selects(self):
        self.assertTrue(
            all(not row["selects_9_over_5"] for row in self.d["candidate_rows"])
        )

    def test_many_rate_time_pairs(self):
        self.assertGreater(self.d["number_of_displayed_rate_time_pairs"], 500)

    def test_no_eta_source(self):
        self.assertFalse(self.d["target_independent_eta_source_exported"])


class TestProgram206(unittest.TestCase):
    def setUp(self):
        self.d = P["206"]

    def test_contract_hash(self):
        path = ROOT / self.d["constructed_object"]["path"]
        self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), self.d["constructed_object"]["sha256"])

    def test_contract_canonical(self):
        verify_canonical_record(ROOT / self.d["constructed_object"]["path"])

    def test_arb_absent(self):
        self.assertFalse(self.d["python_flint_available"])
        self.assertFalse(self.d["arb_binary_available"])

    def test_docker_server_not_operational(self):
        self.assertTrue(self.d["docker_cli_available"])
        self.assertFalse(self.d["docker_server_accessible"])

    def test_no_false_certificate(self):
        self.assertFalse(self.d["formal_five_node_run_executed"])
        self.assertFalse(self.d["formal_sub_three_percent_certificate"])


class TestProgram207(unittest.TestCase):
    def setUp(self):
        self.d = P["207"]

    def test_contract_hash(self):
        path = ROOT / self.d["constructed_object"]["path"]
        self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), self.d["constructed_object"]["sha256"])

    def test_contract_canonical(self):
        verify_canonical_record(ROOT / self.d["constructed_object"]["path"])

    def test_wider_exact_pack(self):
        self.assertEqual(self.d["new_exact_weighted_circulant_cases"], 200)

    def test_source_hash(self):
        path = ROOT / self.d["lean_source"]
        self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), self.d["lean_source_sha256"])
        core = ROOT / self.d["core_probe_source"]
        self.assertEqual(
            hashlib.sha256(core.read_bytes()).hexdigest(),
            self.d["core_probe_source_sha256"],
        )

    def test_scoped_compile_status(self):
        self.assertTrue(self.d["installed_toolchains"])
        self.assertTrue(self.d["core_probe_machine_compiled"])
        self.assertFalse(self.d["machine_compiled"])
        self.assertIn("Mathlib", self.d["compilation_output_tail"])


class TestProgram208(unittest.TestCase):
    def setUp(self):
        self.d = P["208"]

    def test_trace_preserving(self):
        self.assertLess(
            self.d["perfect_catalyst_channel_trace_preservation_residual"], 1e-14
        )

    def test_covariant(self):
        self.assertLess(self.d["perfect_catalyst_covariance_residual"], 1e-14)

    def test_exact_mapping(self):
        self.assertLess(self.d["perfect_catalyst_mapping_residual"], 1e-14)

    def test_partial_scan_fails_necessary_condition(self):
        self.assertLess(self.d["best_partial_catalyst_necessary_margin"], -0.5)

    def test_perfect_saturates(self):
        self.assertLess(abs(self.d["perfect_catalyst_necessary_margin"]), 1e-8)

    def test_arbitrary_target_unlock(self):
        self.assertTrue(self.d["arbitrary_target_conversion_with_perfect_reference"])

    def test_no_selector_source(self):
        self.assertFalse(self.d["strict_selector_source_exported"])


class TestProgram209(unittest.TestCase):
    def setUp(self):
        self.d = P["209"]

    def test_hidden_thinning_nontrivial(self):
        self.assertGreater(self.d["selected_lag"], 100)
        self.assertLess(self.d["thinned_count"], 100)

    def test_coupling_budget(self):
        self.assertLessEqual(self.d["per_sequence_coupling_budget"], 0.05 / 4)

    def test_valid_epsilon_larger(self):
        self.assertGreater(self.d["valid_epsilon"], self.d["nominal_epsilon"])

    def test_nominal_undercoverage(self):
        self.assertLess(self.d["nominal_coverage"], 0.95)

    def test_valid_coverage(self):
        self.assertGreaterEqual(self.d["coupling_valid_coverage"], 0.95)

    def test_information_cost(self):
        self.assertGreater(self.d["valid_mean_width"], self.d["nominal_mean_width"])


class TestProgram210(unittest.TestCase):
    def setUp(self):
        self.d = P["210"]
        self.rows = {row["model"]: row for row in self.d["rows"]}

    def test_ville_threshold(self):
        self.assertAlmostEqual(self.d["threshold"], 20)

    def test_fresh_calibration(self):
        self.assertEqual(self.d["calibration_records_per_step"], 15)

    def test_null_false_alarm(self):
        self.assertLessEqual(
            self.rows["iid_gaussian"]["time_uniform_rejection_rate"], 0.05
        )

    def test_alternative_power(self):
        for model in ["sorted", "AR1", "repeated", "drift"]:
            self.assertGreaterEqual(
                self.rows[model]["time_uniform_rejection_rate"], 0.95
            )

    def test_fast_crossing(self):
        self.assertEqual(self.rows["AR1"]["mean_crossing_time"], 2.0)


class TestProgram211(unittest.TestCase):
    def setUp(self):
        self.d = P["211"]

    def test_regions_nonempty(self):
        self.assertEqual(self.d["nonempty_regions"], self.d["replicates"])

    def test_coverage(self):
        self.assertGreaterEqual(self.d["simultaneous_coverage"], 0.95)

    def test_memory_power(self):
        self.assertGreaterEqual(self.d["memory_detection_power"], 0.95)

    def test_widths_positive(self):
        self.assertTrue(
            all(x > 0 for x in self.d["mean_parameter_interval_widths"].values())
        )

    def test_truth_physical(self):
        truth = self.d["truth"]
        self.assertLess(abs(truth["covariance"]), truth["variance"])
        self.assertLessEqual(truth["blur"], 1)


class TestProgram212(unittest.TestCase):
    def setUp(self):
        self.d = P["212"]

    def test_exact_bounds(self):
        self.assertEqual(self.d["sharp_c2_lower_exact"], "7/25")
        self.assertEqual(self.d["sharp_c2_upper_exact"], "1")

    def test_numeric_bounds(self):
        self.assertAlmostEqual(self.d["sharp_c2_lower"], 0.28)
        self.assertAlmostEqual(self.d["sharp_c2_upper"], 1.0)

    def test_lower_witness(self):
        self.assertAlmostEqual(
            self.d["lower_extremizer_value"], self.d["sharp_c2_lower"]
        )

    def test_upper_witness(self):
        self.assertAlmostEqual(
            self.d["upper_extremizer_value"], self.d["sharp_c2_upper"]
        )

    def test_psd_grid(self):
        self.assertGreater(
            self.d["minimum_psd_eigenvalue_on_interval_grid"], -1e-12
        )


class TestProgram213(unittest.TestCase):
    def setUp(self):
        self.d = P["213"]

    def test_certificate(self):
        path = ROOT / self.d["constructed_object"]["path"]
        self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), self.d["constructed_object"]["sha256"])
        verify_canonical_record(path)

    def test_eight_images(self):
        self.assertEqual(self.d["distinct_natural_endomorphism_images"], 8)

    def test_nogo_hierarchy(self):
        self.assertTrue(self.d["continuous_no_go"])
        self.assertTrue(self.d["borel_homomorphism_no_go"])
        self.assertTrue(self.d["trivial_action_cocycle_no_go"])

    def test_no_phase_source(self):
        self.assertFalse(self.d["strict_phase_source_exported"])


class TestProgram214(unittest.TestCase):
    def setUp(self):
        self.d = P["214"]
        self.rows = {row["challenge"]: row for row in self.d["challenge_rows"]}

    def test_protocol_canonical(self):
        verify_canonical_record(ROOT / self.d["constructed_object"]["path"])

    def test_strict_row_sum(self):
        self.assertAlmostEqual(self.d["strict_row_sum_s"], 1.660307278766099)

    def test_all_scales_pass(self):
        scale_rows = [
            row for row in self.d["challenge_rows"] if row["challenge"].startswith("strict_scale_")
        ]
        self.assertTrue(all(row["accepted"] for row in scale_rows))

    def test_permutation_passes(self):
        self.assertTrue(self.rows["strict_node_permutation"]["accepted"])

    def test_small_noise_passes(self):
        self.assertTrue(self.rows["strict_edge_noise_0.01"]["accepted"])

    def test_large_alternatives_fail(self):
        self.assertFalse(self.rows["strict_edge_noise_0.10"]["accepted"])
        self.assertFalse(self.rows["nearest_neighbour_C12"]["accepted"])

    def test_legacy_fails_structural_gate(self):
        row = self.rows["legacy_C12_raw_generator"]
        self.assertFalse(row["accepted"])
        self.assertFalse(row["structural_gate"]["psd"])


class TestProgram215(unittest.TestCase):
    def setUp(self):
        self.d = P["215"]

    def test_request_canonical(self):
        verify_canonical_record(ROOT / self.d["constructed_object"]["path"])

    def test_no_admitted_bundle(self):
        self.assertEqual(self.d["admitted_bundles"], [])
        self.assertFalse(self.d["external_11_of_11_bundle_found"])

    def test_synthetic_not_external(self):
        row = max(self.d["candidate_rows"], key=lambda x: x["passed_fields"])
        self.assertEqual(row["passed_fields"], 10)
        self.assertEqual(row["source_class"], "synthetic_method_validation")
        self.assertFalse(row["admitted"])

    def test_no_web(self):
        self.assertFalse(self.d["web_or_firecrawl_used"])


class TestProgram216(unittest.TestCase):
    def setUp(self):
        self.d = P["216"]

    def test_lock_canonical(self):
        verify_canonical_record(ROOT / self.d["constructed_object"]["path"])

    def test_gate_failed(self):
        self.assertFalse(self.d["program215_gate_passed"])

    def test_not_executed(self):
        self.assertFalse(self.d["external_prediction_executed"])

    def test_no_validation(self):
        self.assertFalse(self.d["external_physical_validation_claimed"])


if __name__ == "__main__":
    unittest.main()
