#!/usr/bin/env python3
"""Acceptance tests for FIN Release 10.58, programs ST46--ST57."""

from __future__ import annotations

import hashlib
import json
import unittest
from fractions import Fraction
from pathlib import Path

import numpy as np

import fin_st46_st57_research as research


ROOT = Path(__file__).resolve().parent
RESULTS = json.loads((ROOT / "FIN_ST46_ST57_Results.json").read_text(encoding="utf-8"))


def packet_digest(path: Path) -> tuple[str, str]:
    packet = json.loads(path.read_text(encoding="utf-8"))
    key = "configuration" if "configuration" in packet else "specification"
    raw = json.dumps(packet[key], sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(raw).hexdigest(), packet["sha256"]


class TestST46ST57(unittest.TestCase):
    def test_01_inventory(self) -> None:
        self.assertEqual(RESULTS["programs"], "ST46-ST57")
        self.assertTrue(all(f"ST{i}" in RESULTS for i in range(46, 58)))

    def test_02_st46_partial_upstream_certificate(self) -> None:
        result = RESULTS["ST46"]
        lo, hi = result["collision_speed_interval"]
        self.assertLess(lo, 1.278014751820449)
        self.assertGreater(hi, 1.278014751820449)
        self.assertLess(hi - lo, 1e-8)
        self.assertTrue(any(row["strict_inclusion"] for row in result["radius_rows"]))
        self.assertIn("not_full_interval", result["status"])

    def test_03_st46_certificate_matches_results(self) -> None:
        certificate = json.loads((ROOT / "FIN_ST46_Upstream_Sensitivity_Certificate.json").read_text(encoding="utf-8"))
        self.assertEqual(certificate["collision_speed_interval"], RESULTS["ST46"]["collision_speed_interval"])
        self.assertGreater(certificate["memory_operator_error_bound"], certificate["strict_matrix_operator_error_bound"])

    def test_04_st47_trace_law(self) -> None:
        result = RESULTS["ST47"]
        self.assertEqual(result["unique_nonnegative_q_under_trace_density_conservation"], 0.0)
        for row in result["rows"]:
            self.assertAlmostEqual(row["trace_density_minus_coarse"], row["q"], places=12)
            self.assertLess(row["coarse_intertwining_residual"], 1e-12)

    def test_05_st48_gain_sign_is_not_identified(self) -> None:
        result = RESULTS["ST48"]
        self.assertLess(result["reflection_odd_tangent_residual"], 1e-12)
        rows = {row["mu"]: row for row in result["rows"]}
        self.assertLess(rows[-0.35]["final_amplitude"], 1e-8)
        self.assertGreater(rows[0.35]["final_amplitude"], 0.5)
        self.assertTrue(all(row["strict_fixed_point_stationarity"] for row in rows.values()))

    def test_06_st49_projective_joint_object(self) -> None:
        result = RESULTS["ST49"]
        holonomy = complex(result["projective_holonomy"]["real"], result["projective_holonomy"]["imag"])
        self.assertLess(abs(holonomy + 1), 1e-12)
        self.assertEqual(result["density_algebra_dimension"], 144)
        self.assertEqual(result["joint_commutant_dimension"], 1)
        self.assertEqual(result["density_orbit_size"], 24)
        self.assertIn("QW-2191 remains open", result["boundary"])

    def test_07_st49_live_holonomy(self) -> None:
        angles = 2.0 * np.pi * np.arange(research.N) / research.N
        rays = np.stack([np.cos(angles / 2.0), np.sin(angles / 2.0)], axis=1).astype(complex)
        self.assertLess(abs(research.projective_holonomy(rays) + 1), 1e-12)

    def test_08_st50_intertwiner_classification(self) -> None:
        result = RESULTS["ST50"]
        matrix = np.asarray(result["intertwiner_dimension_matrix"])
        self.assertTrue(np.all(np.diag(matrix) == 22))
        self.assertEqual(matrix[0, 1], 1)
        self.assertEqual(matrix[0, 4], 0)
        self.assertEqual(result["invertible_cross_channel_similarity_pairs"], [])
        self.assertTrue(all(result["injective_on_seven_distinct_strict_eigenvalues"].values()))

    def test_09_st51_preregistration_and_power(self) -> None:
        computed, embedded = packet_digest(ROOT / "FIN_ST51_Two_Carrier_Preregistration.json")
        self.assertEqual(computed, embedded)
        self.assertEqual(embedded, RESULTS["ST51"]["preregistration_sha256"])
        self.assertLess(RESULTS["ST51"]["transported_probability_residual"], 2e-12)
        self.assertLess(RESULTS["ST51"]["holdout_false_rejection_rate"], 0.02)
        self.assertEqual(RESULTS["ST51"]["mismatched_detection_power"], 1.0)

    def test_10_st52_exact_root_brackets(self) -> None:
        intervals = RESULTS["ST52"]["critical_y_rational_intervals"]
        for endpoints in intervals.values():
            left, right = map(Fraction, endpoints)
            self.assertLess(left, right)
            self.assertLess(research.polynomial_sign(left) * research.polynomial_sign(right), 0)
        self.assertGreater(RESULTS["ST52"]["radial_hessian_curvature"], 0.4)
        self.assertGreater(RESULTS["ST52"]["phase_hessian_positive_gap"], 0.75)

    def test_11_st53_entropy_ledger(self) -> None:
        result = RESULTS["ST53"]
        self.assertGreater(result["transition_minimum_entry"], 0)
        self.assertLess(result["row_sum_residual"], 1e-12)
        self.assertGreaterEqual(result["shannon_after_strict_markov"], result["shannon_before"])
        self.assertLessEqual(result["relative_entropy_after"], result["relative_entropy_before"])
        self.assertLess(result["unitary_von_neumann_entropy_residual"], 1e-12)
        self.assertGreater(result["normalized_heat_filter_affinity_defect"], 0.1)

    def test_12_st54_universal_property_boundary(self) -> None:
        result = RESULTS["ST54"]
        self.assertIn("unique map", result["universal_property"])
        self.assertIn("classification theorem", result["boundary"])

    def test_13_st55_specification_hash(self) -> None:
        computed, embedded = packet_digest(ROOT / "FIN_ST55_Two_Carrier_Executable_Spec.json")
        self.assertEqual(computed, embedded)
        self.assertEqual(embedded, RESULTS["ST55"]["specification_sha256"])
        self.assertEqual(RESULTS["ST55"]["required_raw_fields"], 7)
        self.assertEqual(RESULTS["ST55"]["logical_role_count"], 4)

    def test_14_st56_compression_nonselection(self) -> None:
        result = RESULTS["ST56"]
        self.assertGreater(result["effective_exponent_spread"], 1.7)
        self.assertTrue(result["directed_ratio_intervals_pairwise_disjoint"])
        self.assertTrue(all(row["coarse_reconstruction_loss"] < 1e-25 for row in result["refinement_rows"]))
        self.assertEqual(result["refinement_rows"][0]["trace_complexity_penalty"], 0.0)

    def test_15_st57_axiom_reconciliation(self) -> None:
        result = RESULTS["ST57"]
        self.assertEqual(result["axiom_count_before"], 9)
        self.assertEqual(result["axiom_count_after_strict_source_audit"], 9)
        self.assertEqual(len(result["rows"]), 9)

    def test_16_recommended_next_batch(self) -> None:
        next_programs = RESULTS["recommended_next_programs"]
        self.assertEqual([row["id"] for row in next_programs], [f"ST{i}" for i in range(58, 70)])

    def test_17_global_boundary(self) -> None:
        boundary = RESULTS["epistemic_boundary"]
        for phrase in ("strict feedback-gain source", "physical scale", "ToE closure"):
            self.assertIn(phrase, boundary)


if __name__ == "__main__":
    unittest.main(verbosity=2)
