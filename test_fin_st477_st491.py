#!/usr/bin/env python3
"""Independent checks for FIN ST477--ST491."""

from __future__ import annotations

import hashlib
import json
import math
import subprocess
import tempfile
import unittest
from pathlib import Path

import numpy as np

import fin_st477_st491_research as research


ROOT = Path(__file__).resolve().parent


class ST477ST491Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r = json.loads((ROOT / "FIN_ST477_ST491_Results.json").read_text())

    def test_program_range_and_packet_hashes(self):
        self.assertEqual(set(self.r), {f"ST{k}" for k in range(477, 492)})
        for k in range(477, 492):
            row = self.r[f"ST{k}"]; path = ROOT / row["packet_file"]
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), row["packet_sha256"])

    def test_exact_two_prime_modular_ranks(self):
        self.assertEqual(self.r["ST477"]["rank"], 1791)
        self.assertEqual(self.r["ST478"]["rank"], 1791)
        self.assertEqual(self.r["ST477"]["nullity"], 237)
        self.assertEqual(self.r["ST478"]["nullity"], 237)
        self.assertEqual(self.r["ST479"]["characteristic_zero_decomposable_rank_interval"], [1791, 1892])
        self.assertFalse(self.r["ST479"]["rank_1791_equality_promoted"])

    def test_live_c_rank_helper_on_small_matrix(self):
        binary = research.compile_modrank()
        prime = 1000003
        matrix = np.eye(7, 11, dtype=np.uint32)
        with tempfile.NamedTemporaryFile(suffix=".bin") as f:
            matrix.tofile(f.name)
            p = subprocess.run([str(binary), f.name, "7", "11", str(prime)], capture_output=True, text=True, timeout=10)
        self.assertEqual(p.returncode, 0)
        self.assertIn("rank=7", p.stdout)

    def test_transition_branch_stabilizer(self):
        row = self.r["ST480"]
        self.assertEqual(row["stabilizer_order"], 2)
        self.assertEqual(row["orbit_size"], 12)
        for endpoint in row["endpoint_checks"]:
            self.assertEqual(len(endpoint["stabilizer"]), 2)
            self.assertGreater(endpoint["unique_max_gap"], 0.89)

    def test_fluctuation_source_accounting(self):
        row = self.r["ST481"]
        self.assertGreater(row["exact_expected_q"], 0)
        self.assertLess(row["relative_sampling_error"], 2e-4)
        self.assertIn("added resource", row["strict_source_status"])
        self.assertFalse(self.r["ST482"]["canonical_selector_exported"])

    def test_dual_threshold_live_equation(self):
        row = self.r["ST483"]; x = row["first_positive_dimensionless_threshold"]
        self.assertLess(abs(math.exp(-x)-2*math.cos(x)), 2e-12)
        self.assertAlmostEqual(research.dual_distance_scalar(x), 1.0, places=13)
        self.assertFalse(row["physical_time"])

    def test_operational_conjugacy_and_clock_boundary(self):
        self.assertLess(self.r["ST484"]["maximum_record_residual"], 2e-13)
        self.assertFalse(self.r["ST484"]["spectrum_alone_sufficient"])
        self.assertFalse(self.r["ST485"]["single_raw_A_unit_closes_all_three_channels"])
        self.assertFalse(self.r["ST485"]["absolute_seconds_generated"])

    def test_transition_identity_not_globalized(self):
        row = self.r["ST486"]
        self.assertLess(row["distance_to_ST342_floating_center"], 2e-14)
        self.assertEqual(row["multistarts_in_same_orbit"], 108)
        self.assertFalse(row["global_first_identity"])

    def test_cover_stop_and_ir_repair(self):
        self.assertFalse(self.r["ST487"]["new_global_lower_theorem"])
        row = self.r["ST488"]
        self.assertGreaterEqual(row["certified_stop"], 0.196-1e-14)
        self.assertGreater(row["minimum_Krawczyk_margin"], 0)
        self.assertLess(row["maximum_weighted_contraction"], 1)

    def test_sampled_degree8_audit_is_not_exhaustive(self):
        self.assertEqual(self.r["ST489"]["true_dependencies_found"], 0)
        self.assertFalse(self.r["ST489"]["exhaustive_theorem"])

    def test_gates(self):
        row = self.r["ST490"]
        self.assertFalse(row["new_strict_gain_source"])
        self.assertFalse(row["new_nonpremise_selector"])
        self.assertFalse(row["new_scale_charged_source"])
        self.assertEqual(self.r["ST491"]["independent_laboratory_record"], "absent")


if __name__ == "__main__":
    unittest.main(verbosity=2)
