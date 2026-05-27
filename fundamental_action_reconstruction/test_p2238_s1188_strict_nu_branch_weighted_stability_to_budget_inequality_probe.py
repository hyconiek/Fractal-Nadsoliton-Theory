from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2238S1188StrictNuBranchWeightedStabilityToBudgetInequalityProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2238_s1188_strict_nu_branch_weighted_stability_to_budget_inequality_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2238_s1188_strict_nu_branch_weighted_stability_to_budget_inequality_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2238_s1188_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["weighted_to_budget_inequality_exported"])
        self.assertTrue(g["weighted_slope_positive"])
        self.assertTrue(g["inequality_holds_against_min_budget"])


if __name__ == "__main__":
    unittest.main()
