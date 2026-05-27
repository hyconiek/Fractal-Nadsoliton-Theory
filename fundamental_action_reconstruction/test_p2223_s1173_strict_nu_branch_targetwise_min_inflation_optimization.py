from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2223S1173StrictNuBranchTargetwiseMinInflationOptimization(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2223_s1173_strict_nu_branch_targetwise_min_inflation_optimization.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2223_s1173_strict_nu_branch_targetwise_min_inflation_optimization.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2223_s1173_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["targetwise_optimization_exported"])
        self.assertTrue(g["all_targetwise_factors_ge_1"])
        self.assertTrue(g["mean_cost_not_worse_than_global"])


if __name__ == "__main__":
    unittest.main()
