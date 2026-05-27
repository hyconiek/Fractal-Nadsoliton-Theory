from __future__ import annotations
import json, subprocess, sys, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2245S1195StrictNuBranchGroupPolicyReserveFactorCalibrationProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2245_s1195_strict_nu_branch_group_policy_reserve_factor_calibration_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2245_s1195_strict_nu_branch_group_policy_reserve_factor_calibration_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2245_s1195_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["reserve_factor_exported"])
        self.assertTrue(g["kappa_capacity_computable"])
        self.assertTrue(g["target_risk_applied"])


if __name__ == "__main__":
    unittest.main()
