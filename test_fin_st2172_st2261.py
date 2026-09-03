#!/usr/bin/env python3
import hashlib
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parent
ROUNDS = [(2172, 2186), (2187, 2201), (2202, 2216), (2217, 2231), (2232, 2246), (2247, 2261)]


class TestST2172ST2261(unittest.TestCase):
    def test_all_ninety_packets_and_hashes(self):
        for lo, hi in ROUNDS:
            results = json.loads((ROOT / f"FIN_ST{lo}_ST{hi}_Results.json").read_text())
            self.assertEqual(list(results), [f"ST{k}" for k in range(lo, hi + 1)])
            self.assertEqual(len(results), 15)
            for k in range(lo, hi + 1):
                row = results[f"ST{k}"]
                digest = hashlib.sha256((ROOT / row["packet_file"]).read_bytes()).hexdigest()
                self.assertEqual(digest, row["packet_sha256"])

    def test_round_one_cycle_and_algebra(self):
        r = json.loads((ROOT / "FIN_ST2172_ST2186_Results.json").read_text())
        self.assertEqual(r["ST2172"]["max_D12_conjugation_residual_inf"], 0.0)
        self.assertEqual(r["ST2180"]["positive"], 220)
        self.assertEqual(r["ST2181"]["orbit_count"], 12)
        self.assertFalse(r["ST2183"]["irreducible_connected_third_order_from_pair_kernel"])

    def test_state_and_adaptive_no_go(self):
        r2 = json.loads((ROOT / "FIN_ST2187_ST2201_Results.json").read_text())
        self.assertAlmostEqual(r2["ST2189"]["energy_plus"], r2["ST2189"]["energy_minus"], places=14)
        self.assertAlmostEqual(r2["ST2190"]["current_plus"], -r2["ST2190"]["current_minus"], places=14)
        self.assertEqual(r2["ST2191"]["mixture_current"], 0.0)
        r3 = json.loads((ROOT / "FIN_ST2202_ST2216_Results.json").read_text())
        self.assertLess(r3["ST2209"]["W_lambda_min"], 0.0)
        self.assertFalse(r3["ST2215"]["adaptive_source_gate_pass"])

    def test_hidden_synergy_exact_counterexample(self):
        r = json.loads((ROOT / "FIN_ST2217_ST2231_Results.json").read_text())
        self.assertEqual(r["ST2220"]["maximum_pair_marginal_difference"], 0.0)
        self.assertEqual(r["ST2221"]["sample_triples"]["0.75"], 0.75)
        self.assertFalse(r["ST2223"]["pair_to_triple_reconstruction_exists"])
        self.assertEqual(r["ST2225"]["coarse_total_variation"], 0.0)

    def test_hodge_counterfamily_and_final_gate(self):
        r5 = json.loads((ROOT / "FIN_ST2232_ST2246_Results.json").read_text())
        self.assertGreater(r5["ST2239"]["eigenvalue_vector_residual_inf"], 20.0)
        self.assertFalse(r5["ST2245"]["unique_strict_model"])
        r6 = json.loads((ROOT / "FIN_ST2247_ST2261_Results.json").read_text())
        self.assertEqual(r6["ST2256"]["strict_mathematical_passes"], 3)
        self.assertEqual(r6["ST2256"]["physical_source_passes"], 0)
        self.assertFalse(r6["ST2261"]["strict_ToE_closure"])
        self.assertTrue((ROOT / r6["ST2261"]["figure"]).exists())


if __name__ == "__main__":
    unittest.main()
