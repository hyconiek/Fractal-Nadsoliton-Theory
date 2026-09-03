#!/usr/bin/env python3
"""Independent checks for the FIN ST462--ST476 packet."""

from __future__ import annotations

import hashlib
import json
import math
import unittest
from pathlib import Path

import numpy as np

import fin_st462_st476_research as research


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "FIN_ST462_ST476_Results.json"


class ST462ST476Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r = json.loads(RESULTS.read_text(encoding="utf-8"))

    def test_complete_program_range_and_packet_hashes(self):
        self.assertEqual(set(self.r), {f"ST{k}" for k in range(462, 477)})
        for k in range(462, 477):
            row = self.r[f"ST{k}"]
            path = ROOT / row["packet_file"]
            self.assertTrue(path.is_file())
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), row["packet_sha256"])

    def test_transition_orbit_is_twelve_fold_after_reflection_quotient(self):
        row = self.r["ST462"]
        self.assertEqual(row["D12_orbit_size"], 12)
        self.assertLess(row["minimum_reflection_stabilizer_residual"], row["orbit_identification_tolerance"])
        self.assertEqual(row["best_orbit_hit_count"], row["start_count"])
        self.assertTrue(2.902496471747767 <= row["best_ratio"] <= 2.9024964917477667)
        self.assertLess(row["maximum_orbit_matching_residual"], 1e-7)

    def test_adversarial_search_did_not_beat_candidate(self):
        row = self.r["ST463"]
        self.assertEqual(row["all_proper_uniform_supports_tested"], 4094)
        self.assertEqual(row["random_Dirichlet_samples"], 60000)
        self.assertFalse(row["candidate_beaten"])
        self.assertGreater(row["minimum_refined_boundary_ratio"], self.r["ST462"]["best_ratio"])

    def test_global_cover_stop_is_not_promoted(self):
        row = self.r["ST464"]
        self.assertFalse(row["cover"]["closed"])
        self.assertFalse(row["new_global_lower_theorem"])
        self.assertLessEqual(row["cover"]["failed_lower_bound"], 0)

    def test_autocatalytic_receiver_needs_a_seed(self):
        row = self.r["ST465"]
        self.assertEqual(row["symmetric_vacuum_simulation"]["maximum_g"], 0)
        self.assertGreater(row["perturbed_simulation"]["maximum_g"], 0)
        self.assertFalse(row["new_autonomous_gain_source_from_symmetric_vacuum"])

    def test_exact_sampled_circuit_audit(self):
        row = self.r["ST466"]
        self.assertEqual(row["true_sampled_circuits"], 0)
        self.assertEqual(sum(x["sampled_subsets"] for a in row["exact_modular_random_audits"] for x in a["rows"]), 72000)

    def test_degree8_exact_attempt_not_overclaimed(self):
        self.assertFalse(self.r["ST467"]["exact_rank_theorem"])
        self.assertIsNone(self.r["ST468"]["exact_rank_returned"])
        self.assertFalse(self.r["ST468"]["explicit_nullspace_returned"])

    def test_boundary_exclusion_and_ir_extension(self):
        self.assertIn("epsilon*log(epsilon)", self.r["ST469"]["dominant_term"])
        row = self.r["ST470"]
        self.assertGreater(row["certified_stop"], 0.1)
        self.assertGreater(row["minimum_Krawczyk_margin"], 0)
        self.assertFalse(row["extension_to_target_0_2"])

    def test_dual_discriminant_live_stationarity(self):
        row = self.r["ST471"]
        x = row["first_positive_global_maximizer"]
        residual = math.sin(x) + math.cos(x) - math.exp(-x)
        self.assertLess(abs(residual), 2e-12)
        self.assertAlmostEqual(research.dual_distance_scalar(x), row["universal_maximum_distance"], places=14)
        self.assertGreater(row["universal_maximum_distance"], 1)

    def test_dual_intertwining_and_isospectral_counterexample(self):
        self.assertTrue(self.r["ST472"]["all_bounds_pass"])
        row = self.r["ST473"]
        self.assertLess(row["spectrum_max_absolute_residual"], 1e-12)
        self.assertLess(row["constant_mode_residual"], 1e-12)
        self.assertGreater(row["heat_vertex_record_max_difference_t_0_8"], 0.02)
        self.assertGreater(row["unitary_Born_vertex_record_max_difference_t_0_8"], 0.04)

    def test_live_ratio_matches_stored_candidate(self):
        p = np.array(self.r["ST462"]["best_probability"], dtype=float)
        self.assertAlmostEqual(research.probability_ratio(p), self.r["ST462"]["best_ratio"], places=13)
        self.assertAlmostEqual(float(np.sum(p)), 1.0, places=13)
        self.assertGreater(float(np.min(p)), 0)

    def test_source_and_evidence_gates_remain_closed(self):
        self.assertFalse(self.r["ST474"]["new_strict_gain_source"])
        self.assertFalse(self.r["ST475"]["new_nonpremise_selector"])
        self.assertFalse(self.r["ST475"]["new_scale_charged_source"])
        self.assertEqual(self.r["ST476"]["independent_laboratory_record"], "absent")


if __name__ == "__main__":
    unittest.main(verbosity=2)
