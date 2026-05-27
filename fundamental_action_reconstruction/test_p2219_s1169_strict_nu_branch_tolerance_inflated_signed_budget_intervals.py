from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2219S1169StrictNuBranchToleranceInflatedSignedBudgetIntervals(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2219_s1169_strict_nu_branch_tolerance_inflated_signed_budget_intervals.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2219_s1169_strict_nu_branch_tolerance_inflated_signed_budget_intervals.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2219_s1169_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["inflated_interval_map_exported"])
        self.assertTrue(g["interval_endpoints_ordered"])
        self.assertTrue(g["above_intervals_include_positive_side"])
        self.assertTrue(g["below_intervals_include_negative_side"])


if __name__ == "__main__":
    unittest.main()
