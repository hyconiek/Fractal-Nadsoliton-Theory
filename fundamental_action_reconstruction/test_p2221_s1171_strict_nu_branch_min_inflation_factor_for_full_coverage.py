from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2221S1171StrictNuBranchMinInflationFactorForFullCoverage(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2221_s1171_strict_nu_branch_min_inflation_factor_for_full_coverage.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2221_s1171_strict_nu_branch_min_inflation_factor_for_full_coverage.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2221_s1171_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["min_inflation_factor_exported"])
        self.assertTrue(g["inflation_factor_ge_1"])
        self.assertTrue(g["all_targets_covered_with_min_factor"])


if __name__ == "__main__":
    unittest.main()
