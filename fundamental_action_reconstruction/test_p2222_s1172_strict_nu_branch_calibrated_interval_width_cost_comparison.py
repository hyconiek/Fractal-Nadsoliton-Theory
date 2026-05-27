from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2222S1172StrictNuBranchCalibratedIntervalWidthCostComparison(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2222_s1172_strict_nu_branch_calibrated_interval_width_cost_comparison.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2222_s1172_strict_nu_branch_calibrated_interval_width_cost_comparison.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2222_s1172_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["width_cost_comparison_exported"])
        self.assertTrue(g["all_width_ratios_ge_1"])
        self.assertTrue(g["mean_width_ratio_ge_1"])


if __name__ == "__main__":
    unittest.main()
