#!/usr/bin/env python3
"""Integrity and semantic regression tests for ST1002--ST1091."""
import hashlib
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parent
ROUNDS = [(1002,1016),(1017,1031),(1032,1046),(1047,1061),(1062,1076),(1077,1091)]


class TestTotalNadsolitonCycle(unittest.TestCase):
    def test_all_packets_and_hashes(self):
        for lo, hi in ROUNDS:
            data = json.loads((ROOT/f"FIN_ST{lo}_ST{hi}_Results.json").read_text())
            self.assertEqual(list(data), [f"ST{k}" for k in range(lo,hi+1)])
            self.assertEqual(len(data), 15)
            for k in range(lo,hi+1):
                packet = data[f"ST{k}"]
                raw = (ROOT/packet["packet_file"]).read_bytes()
                self.assertEqual(hashlib.sha256(raw).hexdigest(), packet["packet_sha256"])

    def test_ontology_correction(self):
        d = json.loads((ROOT/"FIN_ST1002_ST1016_Results.json").read_text())
        self.assertIn("entire fundamental being", d["ST1002"]["types"]["total_nadsoliton"])
        self.assertFalse(d["ST1004"]["zero_representation_exists"])
        self.assertTrue(d["ST1008"]["mass_t3"] > 0.999999999)
        self.assertTrue(d["ST1008"]["pattern_reduction_factor"] < 1.0)

    def test_countermodels_and_no_go_scope(self):
        d = json.loads((ROOT/"FIN_ST1062_ST1076_Results.json").read_text())
        self.assertTrue(d["ST1064"]["survival_mass"] < 1.0)
        self.assertEqual(d["ST1066"]["hitting_time"], 2.0)
        self.assertFalse(d["ST1067"]["strict_gap_protected"])
        z = json.loads((ROOT/"FIN_ST1077_ST1091_Results.json").read_text())
        self.assertIn("conditional", z["ST1077"]["status"])
        self.assertFalse(z["ST1083"].get("strict_total_persistence", False))
        self.assertFalse(z["ST1091"]["strict_ToE_closure"])
        self.assertEqual(z["ST1091"]["programs"], 90)

    def test_no_silent_legacy_or_selector_transfer(self):
        z = json.loads((ROOT/"FIN_ST1077_ST1091_Results.json").read_text())
        self.assertFalse(z["ST1086"]["legacy_role_transfer_used"])
        self.assertFalse(z["ST1086"]["QW2191_discharged"])


if __name__ == "__main__":
    unittest.main()
