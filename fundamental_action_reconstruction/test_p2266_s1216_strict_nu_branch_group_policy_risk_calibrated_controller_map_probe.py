from __future__ import annotations
import json, subprocess, sys, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2266S1216StrictNuBranchGroupPolicyRiskCalibratedControllerMapProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2266_s1216_strict_nu_branch_group_policy_risk_calibrated_controller_map_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2266_s1216_strict_nu_branch_group_policy_risk_calibrated_controller_map_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2266_s1216_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["risk_controller_map_exported"])
        self.assertTrue(g["all_rho_bounded"])
        self.assertTrue(g["all_kappa_scale_ge_one"])


if __name__ == "__main__":
    unittest.main()
