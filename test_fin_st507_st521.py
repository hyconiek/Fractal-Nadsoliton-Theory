#!/usr/bin/env python3
"""Independent checks for FIN ST507--ST521."""

import hashlib,json,unittest
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent


class ST507ST521Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.r=json.loads((ROOT/"FIN_ST507_ST521_Results.json").read_text())
    def test_programs_and_hashes(self):
        self.assertEqual(set(self.r),{f"ST{k}" for k in range(507,522)})
        for k in range(507,522):
            x=self.r[f"ST{k}"];p=ROOT/x["packet_file"];self.assertEqual(hashlib.sha256(p.read_bytes()).hexdigest(),x["packet_sha256"])
    def test_five_prime_bases(self):
        x=self.r["ST507"];p=ROOT/x["combined_basis_file"];self.assertEqual(hashlib.sha256(p.read_bytes()).hexdigest(),x["combined_basis_sha256"])
        self.assertTrue(x["all_new_ranks_1791"]);self.assertTrue(x["all_new_replays_zero"])
        z=np.load(p);self.assertEqual(len(z.files),5);self.assertTrue(all(z[k].shape==(2028,237) for k in z.files))
        self.assertEqual(self.r["ST508"]["pivot_schema_count"],1)
    def test_lift_not_promoted(self):
        self.assertEqual(self.r["ST509"]["reconstructed_count"],0)
        self.assertEqual(self.r["ST510"]["available_count"],0)
        self.assertFalse(self.r["ST511"]["characteristic_zero_rank_1791_closed"])
    def test_equivariant_vacuum_no_go(self):
        self.assertEqual(self.r["ST512"]["fixed_subspace_dimension_under_D12"],0)
        self.assertFalse(self.r["ST512"]["spontaneous_departure_from_exact_u"])
        self.assertEqual(len(self.r["ST513"]["resource_taxonomy"]),3)
    def test_wick_and_dephasing_boundary(self):
        self.assertLess(self.r["ST514"]["matrix_residual"],1e-12);self.assertFalse(self.r["ST514"]["operational_equivalence"])
        self.assertEqual(self.r["ST515"]["diagonal_mode_factor_residual"],0)
        self.assertGreater(self.r["ST516"]["L1_difference"],.6)
    def test_dual_tester_margin(self):
        x=self.r["ST517"];self.assertAlmostEqual(x["certified_triangle_lower_bound"],.84);self.assertTrue(x["positive_margin"]);self.assertFalse(x["finite_shot_decision_rule"])
    def test_ir_and_curvature(self):
        self.assertGreater(self.r["ST518"]["certified_stop"],.23);self.assertGreater(self.r["ST518"]["minimum_margin"],0)
        self.assertGreater(self.r["ST519"]["minimum_paid_curvature"],5.6);self.assertFalse(self.r["ST519"]["global_minimum_theorem"])
    def test_gates(self):
        x=self.r["ST520"];self.assertFalse(x["new_strict_gain_source"]);self.assertFalse(x["new_nonpremise_selector"]);self.assertFalse(x["new_scale_charged_source"]);self.assertFalse(x["characteristic_zero_degree8_rank_exact"])
        self.assertEqual(self.r["ST521"]["independent_laboratory_record"],"absent")

if __name__=="__main__":unittest.main(verbosity=2)
