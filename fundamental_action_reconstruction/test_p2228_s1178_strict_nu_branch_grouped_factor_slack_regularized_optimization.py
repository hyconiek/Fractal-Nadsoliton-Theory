from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2228S1178StrictNuBranchGroupedFactorSlackRegularizedOptimization(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2228_s1178_strict_nu_branch_grouped_factor_slack_regularized_optimization.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2228_s1178_strict_nu_branch_grouped_factor_slack_regularized_optimization.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2228_s1178_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["grouped_slack_optimization_exported"])
        self.assertTrue(g["candidates_nonempty"])
        self.assertTrue(g["has_fully_feasible_candidate"])


if __name__ == "__main__":
    unittest.main()
