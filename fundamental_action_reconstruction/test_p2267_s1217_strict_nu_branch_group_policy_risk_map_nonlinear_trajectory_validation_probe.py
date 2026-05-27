from __future__ import annotations
import json, subprocess, sys, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2267S1217StrictNuBranchGroupPolicyRiskMapNonlinearTrajectoryValidationProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2267_s1217_strict_nu_branch_group_policy_risk_map_nonlinear_trajectory_validation_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2267_s1217_strict_nu_branch_group_policy_risk_map_nonlinear_trajectory_validation_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2267_s1217_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["risk_map_validation_exported"])
        self.assertTrue(g["all_pass_rates_bounded"])
        self.assertTrue(g["all_targets_bounded"])


if __name__ == "__main__":
    unittest.main()
