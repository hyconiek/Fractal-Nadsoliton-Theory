from __future__ import annotations
import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2271S1221StrictNuBranchGroupPolicySymbolicVsEmpiricalDerivativeSlackProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2271_s1221_strict_nu_branch_group_policy_symbolic_vs_empirical_derivative_slack_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2271_s1221_strict_nu_branch_group_policy_symbolic_vs_empirical_derivative_slack_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2271_s1221_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["slack_rows_exported"])
        self.assertTrue(g["all_symbolic_rho_bounds_cover_empirical"])
        self.assertTrue(g["all_symbolic_kappa_bounds_cover_empirical"])


if __name__ == "__main__":
    unittest.main()
