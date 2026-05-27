from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2227S1177StrictNuBranchGroupedPolicyDirectCoverageMatrix(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2227_s1177_strict_nu_branch_grouped_policy_direct_coverage_matrix.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2227_s1177_strict_nu_branch_grouped_policy_direct_coverage_matrix.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2227_s1177_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["coverage_matrix_exported"])
        self.assertTrue(g["all_rows_covered"])
        self.assertTrue(g["nonnegative_slack_all_rows"])


if __name__ == "__main__":
    unittest.main()
