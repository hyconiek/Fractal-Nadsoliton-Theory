from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2226S1176StrictNuBranchGroupedPolicyWidthCoverageTradeoff(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2226_s1176_strict_nu_branch_grouped_policy_width_coverage_tradeoff.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2226_s1176_strict_nu_branch_grouped_policy_width_coverage_tradeoff.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2226_s1176_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["grouped_tradeoff_exported"])
        self.assertTrue(g["grouped_factors_ge_1"])
        self.assertTrue(g["grouped_not_worse_than_global"])


if __name__ == "__main__":
    unittest.main()
