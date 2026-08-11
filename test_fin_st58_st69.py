#!/usr/bin/env python3
"""Acceptance tests for FIN ST58--ST69."""

from __future__ import annotations

import hashlib
import json
import math
import unittest
from pathlib import Path

import numpy as np

from fin_st58_st69_research import st67_pi_holonomy_chirality


ROOT = Path(__file__).resolve().parent
RESULTS = json.loads((ROOT / "FIN_ST58_ST69_Results.json").read_text(encoding="utf-8"))


def packet_digest(path: Path) -> tuple[str, str]:
    packet = json.loads(path.read_text(encoding="utf-8"))
    payload = packet["specification"]
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(raw).hexdigest(), packet["sha256"]


def all_finite(value) -> bool:
    if isinstance(value, dict):
        return all(all_finite(item) for item in value.values())
    if isinstance(value, list):
        return all(all_finite(item) for item in value)
    if isinstance(value, float):
        return math.isfinite(value)
    return True


class TestST58ST69(unittest.TestCase):
    def test_01_inventory(self) -> None:
        self.assertEqual([RESULTS[f"ST{i}"]["program"] for i in range(58, 70)], [f"ST{i}" for i in range(58, 70)])

    def test_02_st58_full_interval_collision(self) -> None:
        result = RESULTS["ST58"]
        lo, hi = result["collision_speed_interval"]
        self.assertLess(hi - lo, 1.5e-11)
        self.assertLess(lo, 1.278014751820449)
        self.assertGreater(hi, 1.278014751820449)
        self.assertGreater(result["stationary_root"]["minimum_inclusion_margin"], 0.0)
        self.assertLess(result["stationary_root"]["defect_infinity_norm_upper"], 1e-10)
        self.assertTrue(all(row[0] > 0 for row in result["residue_intervals"]))

    def test_03_st58_certificate_matches_results(self) -> None:
        certificate = json.loads((ROOT / "FIN_ST58_Full_Interval_Certificate.json").read_text(encoding="utf-8"))
        self.assertEqual(certificate["collision_speed_interval"], RESULTS["ST58"]["collision_speed_interval"])
        self.assertEqual(certificate["coefficient_intervals"], RESULTS["ST58"]["coefficient_intervals"])

    def test_04_st59_sign_flip_no_go(self) -> None:
        result = RESULTS["ST59"]
        self.assertEqual({row["odd_gradient_flow_gain"] for row in result["rows"]}, {-1.0, 1.0})
        self.assertTrue(all(row["stationary_at_A"] and row["stabilizer_invariant"] for row in result["rows"]))

    def test_05_st60_projective_source_obstruction(self) -> None:
        result = RESULTS["ST60"]
        self.assertTrue(result["strict_spectral_projector_holonomies_all_trivial_when_defined"])
        self.assertAlmostEqual(result["inserted_half_angle_holonomy"]["real"], -1.0, places=12)
        self.assertEqual(sum(not row["defined"] for row in result["integer_mode_rows"]), 1)

    def test_06_st61_heat_signature_no_go(self) -> None:
        result = RESULTS["ST61"]
        lo, hi = result["normalized_second_moment_difference_at_q0_interval"]
        self.assertLess(hi, 0.0)
        self.assertTrue(result["interval_excludes_zero"])
        self.assertLess(hi - lo, 5e-13)

    def test_07_st62_information_count_bracket(self) -> None:
        result = RESULTS["ST62"]
        self.assertEqual(result["necessary_shots_per_preparation_information_bound"], 2)
        self.assertEqual(result["sufficient_shots_per_preparation_chernoff_bound"], 5)
        self.assertGreater(result["ST51_shots_per_preparation"], result["sufficient_shots_per_preparation_chernoff_bound"])

    def test_08_st62_packet_matches(self) -> None:
        packet = json.loads((ROOT / "FIN_ST62_Finite_Count_Bounds.json").read_text(encoding="utf-8"))
        self.assertEqual(packet["chernoff_information_per_shot_per_preparation"], RESULTS["ST62"]["chernoff_information_per_shot_per_preparation"])

    def test_09_st63_cp_pinching(self) -> None:
        result = RESULTS["ST63"]
        self.assertEqual(result["spectral_block_multiplicities"], [1, 2, 2, 2, 2, 2, 1])
        self.assertEqual(result["pinching_output_algebra_dimension"], 22)
        self.assertEqual(result["pinching_kraus_rank"], 7)
        self.assertGreater(result["choi_minimum_eigenvalue"], -2e-14)
        self.assertFalse(result["invertible_CP_cross_intertwiner_exists"])

    def test_10_st64_thermodynamic_identity(self) -> None:
        result = RESULTS["ST64"]
        self.assertLess(result["free_energy_identity_residual"], 1e-13)
        self.assertLess(result["entropy_production_identity_residual"], 1e-13)
        self.assertEqual(len(result["minimal_dimensional_inputs"]), 2)

    def test_11_st65_status_is_not_promoted(self) -> None:
        result = RESULTS["ST65"]
        self.assertEqual(result["localized_hits_at_kappa_1"], 0)
        self.assertLess(result["last_regular_numerical_continuation_kappa"], 0.02)
        self.assertTrue(result["status"].startswith("strong_numerical"))

    def test_12_st66_polynomial_branch_classification(self) -> None:
        result = RESULTS["ST66"]
        self.assertEqual(result["positive_stable_radial_roots"], 1)
        self.assertEqual(result["positive_angular_saddle_radial_roots"], 1)
        self.assertEqual(result["stable_branch_count"], 12)
        self.assertGreater(result["stable_radial_curvature"], 0.0)
        self.assertGreater(result["stable_angular_curvature"], 0.0)

    def test_13_st67_chirality_separation(self) -> None:
        result = RESULTS["ST67"]
        self.assertAlmostEqual(result["forward_holonomy"]["real"], -1.0, places=12)
        self.assertAlmostEqual(result["reverse_holonomy"]["real"], -1.0, places=12)
        self.assertGreater(result["forward_bargmann_chirality"], 0.9)
        self.assertLess(result["reverse_bargmann_chirality"], -0.9)
        self.assertLess(result["chirality_reflection_residual"], 1e-12)

    def test_14_st67_live_replay(self) -> None:
        live = st67_pi_holonomy_chirality()
        self.assertAlmostEqual(live["forward_bargmann_chirality"], RESULTS["ST67"]["forward_bargmann_chirality"], places=14)

    def test_15_st68_validator_hash_and_rejections(self) -> None:
        computed, embedded = packet_digest(ROOT / "FIN_ST68_Calibration_Custody_Validator.json")
        self.assertEqual(computed, embedded)
        result = RESULTS["ST68"]
        self.assertEqual(embedded, result["specification_sha256"])
        self.assertTrue(result["valid_case_accepted"])
        self.assertTrue(result["duplicate_role_case_rejected"])
        self.assertTrue(result["missing_record_case_rejected"])
        self.assertLess(result["actual_endpoint_TV_maximum"], result["analytic_TV_bound"])

    def test_16_st69_axiom_reconciliation(self) -> None:
        result = RESULTS["ST69"]
        self.assertEqual(result["strict_source_group_count_before"], 9)
        self.assertEqual(result["strict_source_group_count_after"], 9)
        self.assertEqual(len(result["rows"]), 9)

    def test_17_recommended_next_batch(self) -> None:
        next_programs = RESULTS["recommended_next_programs"]
        self.assertEqual(len(next_programs), 12)
        self.assertEqual(next_programs[0]["id"], "ST70")
        self.assertEqual(next_programs[-1]["id"], "ST81")

    def test_18_no_nonfinite_serialized_numbers(self) -> None:
        self.assertTrue(all_finite(RESULTS))

    def test_19_figures_exist(self) -> None:
        figures = list((ROOT / "FIN_ST58_ST69_Figures").glob("*.png"))
        self.assertEqual(len(figures), 7)
        self.assertTrue(all(path.stat().st_size > 10_000 for path in figures))

    def test_20_global_boundary(self) -> None:
        boundary = RESULTS["epistemic_boundary"]
        for phrase in ["QW-2191", "laboratory evidence", "legacy role transfer", "ToE closure"]:
            self.assertIn(phrase, boundary)


if __name__ == "__main__":
    unittest.main(verbosity=2)
