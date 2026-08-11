#!/usr/bin/env python3
"""Live acceptance tests for FIN ST70--ST81."""

from __future__ import annotations

import base64
import hashlib
import json
import math
import unittest
from fractions import Fraction
from pathlib import Path

import numpy as np
from cryptography.hazmat.primitives.asymmetric.ed25519 import Ed25519PublicKey

from fin_st01_st15_research import strict_operator
from fin_st70_st81_research import (
    canonical_bytes,
    canonical_digest,
    st71_time_oriented_response_classification,
    st72_spin_central_extension_obstruction,
    st79_orientation_odd_bargmann_obstruction,
)


ROOT = Path(__file__).resolve().parent
RESULTS = json.loads((ROOT / "FIN_ST70_ST81_Results.json").read_text(encoding="utf-8"))


def all_finite(value) -> bool:
    if isinstance(value, dict):
        return all(all_finite(item) for item in value.values())
    if isinstance(value, list):
        return all(all_finite(item) for item in value)
    if isinstance(value, float):
        return math.isfinite(value)
    return True


class TestST70ST81(unittest.TestCase):
    def test_01_inventory(self) -> None:
        self.assertEqual(
            [RESULTS[f"ST{i}"]["program"] for i in range(70, 82)],
            [f"ST{i}" for i in range(70, 82)],
        )

    def test_02_st70_all_exact_checks_pass(self) -> None:
        result = RESULTS["ST70"]
        self.assertTrue(result["all_checks_pass"])
        self.assertTrue(all(result["checks"].values()))

    def test_03_st70_exact_root_sign_and_monotonicity(self) -> None:
        result = RESULTS["ST70"]
        left = list(map(Fraction, result["left_polynomial_interval"]))
        right = list(map(Fraction, result["right_polynomial_interval"]))
        derivative = list(map(Fraction, result["derivative_interval"]))
        self.assertGreater(left[0], 0)
        self.assertLess(right[1], 0)
        self.assertLess(derivative[1], 0)

    def test_04_st70_packet_matches(self) -> None:
        packet = json.loads((ROOT / "FIN_ST70_Rational_Replay_Certificate.json").read_text(encoding="utf-8"))
        self.assertEqual(packet["collision_speed_interval"], RESULTS["ST70"]["collision_speed_interval"])
        self.assertEqual(packet["checks"], RESULTS["ST70"]["checks"])

    def test_05_st71_time_oriented_gain_classification(self) -> None:
        result = RESULTS["ST71"]
        self.assertEqual(result["reflection_odd_state_dimension"], 5)
        for row in result["rows"]:
            self.assertEqual(row["uniform_positive_gain"], row["odd_sector_curvature"] < 0)

    def test_06_st71_live_replay(self) -> None:
        _, a, _ = strict_operator()
        live = st71_time_oriented_response_classification(a)
        self.assertTrue(np.allclose(live["strict_odd_sector_eigenvalues"], RESULTS["ST71"]["strict_odd_sector_eigenvalues"]))

    def test_07_st72_nontrivial_cocycle(self) -> None:
        result = RESULTS["ST72"]
        self.assertEqual(result["cocycle_defects"], 0)
        self.assertGreater(result["coboundary_augmented_rank"], result["coboundary_matrix_rank"])
        self.assertTrue(result["cocycle_is_nontrivial"])

    def test_08_st72_live_rank_replay(self) -> None:
        _, a, _ = strict_operator()
        live = st72_spin_central_extension_obstruction(a)
        self.assertEqual(live["coboundary_matrix_rank"], 10)
        self.assertEqual(live["coboundary_augmented_rank"], 11)

    def test_09_st73_fine_data_selector(self) -> None:
        result = RESULTS["ST73"]
        self.assertAlmostEqual(result["fine_statistic_slope_in_q"], 2.0, places=12)
        self.assertEqual(len(result["online_trajectory"]), 60)
        self.assertTrue(all(row["rate_distortion_selected_q"] <= row["supplied_fine_q"] + 1e-14 for row in result["rows"]))

    def test_10_st74_robust_count_ordering(self) -> None:
        result = RESULTS["ST74"]
        self.assertGreater(result["robust_joint_L1_separation"], 0)
        self.assertLessEqual(result["necessary_total_shots_for_selected_adversarial_pair"], result["distribution_free_sufficient_total_shots"])
        self.assertEqual(result["necessary_mean_shots_per_preparation_ceiling"], 3)
        self.assertEqual(result["distribution_free_sufficient_mean_shots_per_preparation_ceiling"], 1192)

    def test_11_st74_packet_matches(self) -> None:
        packet = json.loads((ROOT / "FIN_ST74_Robust_Count_Design.json").read_text(encoding="utf-8"))
        self.assertEqual(packet["robust_joint_L1_separation"], RESULTS["ST74"]["robust_joint_L1_separation"])

    def test_12_st75_complete_intertwiner_space(self) -> None:
        result = RESULTS["ST75"]
        self.assertEqual(result["fixed_block_algebra_dimension"], 22)
        self.assertEqual(result["complex_linear_intertwiner_space_dimension"], 484)
        self.assertTrue(result["below_two_pi"])
        self.assertTrue(all(row["maximum_intertwining_residual"] < 1e-12 for row in result["example_channels"]))

    def test_13_st76_removal_resources(self) -> None:
        result = RESULTS["ST76"]
        self.assertEqual(result["resource_count"], 4)
        self.assertEqual(len(result["resources"]), 4)
        self.assertTrue(all(row["removal_countermodel"] for row in result["resources"]))

    def test_14_st77_fold_and_return(self) -> None:
        result = RESULTS["ST77"]
        self.assertAlmostEqual(result["fold_kappa"], 0.0165890523863829, places=13)
        self.assertLess(result["fold_system_residual_inf"], 1e-12)
        self.assertGreater(result["augmented_fold_jacobian_minimum_singular_value"], 0.05)
        self.assertTrue(result["return_branch_found"])

    def test_15_st77_not_promoted_to_interval_theorem(self) -> None:
        result = RESULTS["ST77"]
        self.assertTrue(result["status"].startswith("strong_numerical"))
        self.assertIn("does not interval-certify", result["boundary"])

    def test_16_st78_backreaction_threshold(self) -> None:
        result = RESULTS["ST78"]
        self.assertGreater(result["critical_coupling_mu_over_lambda1"], 0.46)
        self.assertLess(result["critical_coupling_mu_over_lambda1"], 0.47)
        unit_row = result["rows"][-1]
        self.assertEqual(unit_row["coupling"], 1.0)
        self.assertEqual(unit_row["positive_radial_root_count"], 0)

    def test_17_st79_chirality_obstruction(self) -> None:
        result = RESULTS["ST79"]
        self.assertLess(result["maximum_strict_candidate_absolute_chirality"], 1e-40)
        self.assertTrue(all(row["odd_sum_residual"] < 1e-40 for row in result["rows"]))

    def test_18_st79_live_replay(self) -> None:
        _, a, _ = strict_operator()
        live = st79_orientation_odd_bargmann_obstruction(a)
        self.assertAlmostEqual(live["maximum_strict_candidate_absolute_chirality"], RESULTS["ST79"]["maximum_strict_candidate_absolute_chirality"], places=55)

    def test_19_st80_signed_packet(self) -> None:
        output = json.loads((ROOT / "FIN_ST80_Signed_Custody_Validator.json").read_text(encoding="utf-8"))
        packet = output["specification"]
        self.assertEqual(canonical_digest(packet), output["packet_sha256"])
        core_hash = canonical_digest(packet["core"])
        for attestation in packet["attestations"]:
            message = attestation["message"]
            self.assertEqual(message["core_hash"], core_hash)
            public = base64.b64decode(packet["core"]["public_keys"][message["role"]])
            Ed25519PublicKey.from_public_bytes(public).verify(
                base64.b64decode(attestation["signature_base64"]), canonical_bytes(message)
            )

    def test_20_st80_rejection_cases(self) -> None:
        result = RESULTS["ST80"]
        self.assertTrue(result["valid_case_accepted"])
        self.assertTrue(result["tampered_event_rejected"])
        self.assertTrue(result["duplicate_public_key_rejected"])

    def test_21_st81_sources_remain(self) -> None:
        result = RESULTS["ST81"]
        self.assertEqual(result["strict_source_group_count_before"], 9)
        self.assertEqual(result["strict_source_group_count_after"], 9)
        self.assertEqual(len(result["rows"]), 9)

    def test_22_next_batch(self) -> None:
        next_programs = RESULTS["recommended_next_programs"]
        self.assertEqual(len(next_programs), 12)
        self.assertEqual(next_programs[0]["id"], "ST82")
        self.assertEqual(next_programs[-1]["id"], "ST93")

    def test_23_figures_and_finite_serialization(self) -> None:
        figures = list((ROOT / "FIN_ST70_ST81_Figures").glob("*.png"))
        self.assertEqual(len(figures), 8)
        self.assertTrue(all(path.stat().st_size > 10_000 for path in figures))
        self.assertTrue(all_finite(RESULTS))

    def test_24_global_boundary(self) -> None:
        boundary = RESULTS["epistemic_boundary"]
        for phrase in ["QW-2191", "laboratory evidence", "legacy role transfer", "ToE closure"]:
            self.assertIn(phrase, boundary)


if __name__ == "__main__":
    unittest.main(verbosity=2)
