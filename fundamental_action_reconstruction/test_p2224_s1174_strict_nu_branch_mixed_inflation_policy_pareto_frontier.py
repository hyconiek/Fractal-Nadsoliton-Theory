from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2224S1174StrictNuBranchMixedInflationPolicyParetoFrontier(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2224_s1174_strict_nu_branch_mixed_inflation_policy_pareto_frontier.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2224_s1174_strict_nu_branch_mixed_inflation_policy_pareto_frontier.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2224_s1174_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["pareto_frontier_exported"])
        self.assertTrue(g["frontier_nonempty"])
        self.assertTrue(g["feasible_point_exists"])


if __name__ == "__main__":
    unittest.main()
