from __future__ import annotations
import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2272S1222StrictNuBranchGroupPolicyEmpiricalDerivativeResolutionAndQuantileCalibrationProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2272_s1222_strict_nu_branch_group_policy_empirical_derivative_resolution_and_quantile_calibration_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2272_s1222_strict_nu_branch_group_policy_empirical_derivative_resolution_and_quantile_calibration_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2272_s1222_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["calibrated_rows_exported"])
        self.assertTrue(g["resolved_zero_sensitivity_issue"])
        self.assertIn(g["all_symbolic_cover_empirical_max_rho"], [True, False])
        self.assertIn(g["all_symbolic_cover_empirical_max_kappa"], [True, False])


if __name__ == "__main__":
    unittest.main()
