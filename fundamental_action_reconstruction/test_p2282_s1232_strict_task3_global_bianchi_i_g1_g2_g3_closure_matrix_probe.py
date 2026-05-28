from __future__ import annotations
import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2282S1232StrictTask3GlobalBianchiIG1G2G3ClosureMatrixProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2282_s1232_strict_task3_global_bianchi_i_g1_g2_g3_closure_matrix_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2282_s1232_strict_task3_global_bianchi_i_g1_g2_g3_closure_matrix_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2282_s1232_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["gap_rows_exported"])
        self.assertTrue(g["contains_G1"])
        self.assertTrue(g["contains_G2"])
        self.assertTrue(g["contains_G3"])
        self.assertTrue(g["closure_score_bounded"])
        rows = data["strict_task3_global_bianchi_i_g1_g2_g3_closure_matrix_probe"]["gap_rows"]
        by_id = {row["id"]: row for row in rows}
        self.assertLess(by_id["G2_nonlinear_trajectory_realism"]["metric"], 1.0)
        self.assertEqual(by_id["G3_operational_policy_rule"]["status"], "OPEN")


if __name__ == "__main__":
    unittest.main()
