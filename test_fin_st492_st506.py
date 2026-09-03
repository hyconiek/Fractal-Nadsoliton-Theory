#!/usr/bin/env python3
"""Independent checks for FIN ST492--ST506."""

from __future__ import annotations

import hashlib
import json
import math
import unittest
from pathlib import Path

import numpy as np

import fin_st492_st506_research as research


ROOT=Path(__file__).resolve().parent


class ST492ST506Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls): cls.r=json.loads((ROOT/"FIN_ST492_ST506_Results.json").read_text())

    def test_program_range_and_hashes(self):
        self.assertEqual(set(self.r),{f"ST{k}" for k in range(492,507)})
        for k in range(492,507):
            row=self.r[f"ST{k}"];path=ROOT/row["packet_file"]
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(),row["packet_sha256"])

    def test_rank_and_pivot_robustness(self):
        a=self.r["ST492"];b=self.r["ST493"]
        self.assertTrue(a["all_ranks_equal_1791"]);self.assertTrue(a["all_pivot_schemas_equal"])
        self.assertTrue(b["all_ranks_equal_1791"]);self.assertTrue(b["alternative_pivots_equal_across_primes"])
        self.assertTrue(self.r["ST494"]["same_schema_between_point_ensembles"])
        self.assertFalse(self.r["ST494"]["rational_schema_theorem"])

    def test_modular_nullspaces_live_replay(self):
        row=self.r["ST495"];path=ROOT/row["nullspace_file"]
        self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(),row["nullspace_sha256"])
        z=np.load(path)
        for prime in (1000003,1000033):
            X=z[str(prime)].astype(np.int64)
            self.assertEqual(X.shape,(2028,237))
            M=research.degree8_matrix_seed(prime,research.SEED+447,npts=8).astype(np.int64)
            self.assertEqual(int(np.count_nonzero((M@X)%prime)),0)

    def test_two_prime_rational_probe_not_promoted(self):
        row=self.r["ST496"]
        self.assertFalse(row["all_reconstructed"])
        self.assertFalse(row["exact_characteristic_zero_polynomial_checks"])
        self.assertTrue(all(not x["reconstructed"] for x in row["reconstruction_rows"]))

    def test_transition_local_curvature_and_stress(self):
        self.assertGreater(self.r["ST497"]["minimum_eigenvalue"],1e-3)
        self.assertEqual(self.r["ST497"]["Morse_index"],0)
        row=self.r["ST498"]
        self.assertEqual(row["run_count"],160);self.assertEqual(row["best_orbit_value_hits"],160)
        self.assertFalse(row["candidate_beaten"]);self.assertTrue(math.isfinite(row["best_ratio"]))

    def test_ir_beyond_point_two(self):
        row=self.r["ST499"]
        self.assertGreater(row["certified_stop"],.2)
        self.assertGreater(row["minimum_margin"],0);self.assertLess(row["maximum_contraction"],1)

    def test_source_and_selector_accounting(self):
        self.assertFalse(self.r["ST500"]["strict_gain_source"])
        self.assertFalse(self.r["ST501"]["strict_selector"])
        self.assertEqual(self.r["ST501"]["QW_2191"],"open")

    def test_dual_probe_and_multitime_conjugacy(self):
        self.assertAlmostEqual(self.r["ST502"]["response_distance"],1.0,places=13)
        self.assertFalse(self.r["ST502"]["physical_tester"])
        self.assertLess(self.r["ST503"]["maximum_residual"],1e-12)
        self.assertFalse(self.r["ST503"]["fixed_untransported_apparatus_equivalent"])

    def test_relative_clock_bounds(self):
        row=self.r["ST504"]
        self.assertLess(row["unitary_heat_ratio_error"],row["unitary_heat_ratio_bound"])
        self.assertLess(row["wave_sqrt_ratio_error"],row["wave_sqrt_ratio_bound"])
        self.assertFalse(row["absolute_clock"])

    def test_gates(self):
        row=self.r["ST505"]
        self.assertFalse(row["new_strict_gain_source"]);self.assertFalse(row["new_nonpremise_selector"])
        self.assertFalse(row["new_scale_charged_source"]);self.assertFalse(row["characteristic_zero_degree8_rank_exact"])
        self.assertEqual(self.r["ST506"]["independent_laboratory_record"],"absent")


if __name__=="__main__":unittest.main(verbosity=2)
