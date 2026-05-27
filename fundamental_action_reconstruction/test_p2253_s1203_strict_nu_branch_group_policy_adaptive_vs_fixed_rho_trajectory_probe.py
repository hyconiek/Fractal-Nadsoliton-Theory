from __future__ import annotations
import json, subprocess, sys, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2253S1203StrictNuBranchGroupPolicyAdaptiveVsFixedRhoTrajectoryProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2253_s1203_strict_nu_branch_group_policy_adaptive_vs_fixed_rho_trajectory_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2253_s1203_strict_nu_branch_group_policy_adaptive_vs_fixed_rho_trajectory_probe.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(data["schema_version"], "p2253_s1203_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["adaptive_vs_fixed_exported"])
        self.assertTrue(g["fixed_nonnegative_margin"])
        self.assertTrue(g["adaptive_nonnegative_margin"])


if __name__ == "__main__":
    unittest.main()
