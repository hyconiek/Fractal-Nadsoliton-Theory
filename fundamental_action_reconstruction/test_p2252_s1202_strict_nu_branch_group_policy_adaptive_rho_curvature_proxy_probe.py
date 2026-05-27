from __future__ import annotations
import json, subprocess, sys, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2252S1202StrictNuBranchGroupPolicyAdaptiveRhoCurvatureProxyProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2252_s1202_strict_nu_branch_group_policy_adaptive_rho_curvature_proxy_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2252_s1202_strict_nu_branch_group_policy_adaptive_rho_curvature_proxy_probe.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(data["schema_version"], "p2252_s1202_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["adaptive_rho_exported"])
        self.assertTrue(g["curvature_proxy_bounded"])
        self.assertTrue(g["rho_within_bounds"])


if __name__ == "__main__":
    unittest.main()
