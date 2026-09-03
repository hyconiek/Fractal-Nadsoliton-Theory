#!/usr/bin/env python3
import hashlib
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parent
ROUNDS = [(2262,2276),(2277,2291),(2292,2306),(2307,2321),(2322,2336),(2337,2351)]


class TestST2262ST2351(unittest.TestCase):
    def test_ninety_packets(self):
        for lo, hi in ROUNDS:
            d = json.loads((ROOT/f"FIN_ST{lo}_ST{hi}_Results.json").read_text())
            self.assertEqual(list(d), [f"ST{k}" for k in range(lo, hi+1)])
            self.assertEqual(len(d), 15)
            for k in range(lo, hi+1):
                row=d[f"ST{k}"]
                self.assertEqual(hashlib.sha256((ROOT/row["packet_file"]).read_bytes()).hexdigest(), row["packet_sha256"])

    def test_A3_A4_fail(self):
        a=json.loads((ROOT/"FIN_ST2262_ST2276_Results.json").read_text())
        self.assertEqual(a["ST2263"]["parallel_energy_additivity_residual"],0.0)
        self.assertFalse(a["ST2274"]["separate_linearity_strict"])
        b=json.loads((ROOT/"FIN_ST2277_ST2291_Results.json").read_text())
        self.assertEqual(b["ST2280"]["Aut_W"],"D12")
        self.assertEqual(b["ST2284"]["equivariant_face_scalar_dimension"],12)
        self.assertFalse(b["ST2289"]["universality_strict"])

    def test_power_theorems(self):
        a=json.loads((ROOT/"FIN_ST2292_ST2306_Results.json").read_text())
        self.assertLess(a["ST2295"]["residual_inf"],1e-12)
        self.assertEqual(a["ST2304"]["p_selected"],1)
        self.assertFalse(a["ST2304"]["strict_selection"])
        b=json.loads((ROOT/"FIN_ST2307_ST2321_Results.json").read_text())
        self.assertEqual(b["ST2309"]["self_similarity_up_to_orbit_permutation_requires_p"],1)
        self.assertFalse(b["ST2313"]["nontrivial_exact_power_self_similarity"])

    def test_instrument_and_gate(self):
        a=json.loads((ROOT/"FIN_ST2322_ST2336_Results.json").read_text())
        self.assertEqual(a["ST2324"]["pair_only_Fisher_information"],0)
        self.assertEqual(a["ST2328"]["sufficient_events"],19)
        self.assertFalse(a["ST2335"]["strict_source_bridge"])
        b=json.loads((ROOT/"FIN_ST2337_ST2351_Results.json").read_text())
        self.assertEqual(b["ST2344"]["mathematical_pass"],7)
        self.assertEqual(b["ST2344"]["strict_source_pass"],3)
        self.assertEqual(b["ST2344"]["physical_pass"],0)
        self.assertFalse(b["ST2351"]["strict_ToE_closure"])
        self.assertTrue((ROOT/b["ST2351"]["figure"]).exists())


if __name__=="__main__": unittest.main()
