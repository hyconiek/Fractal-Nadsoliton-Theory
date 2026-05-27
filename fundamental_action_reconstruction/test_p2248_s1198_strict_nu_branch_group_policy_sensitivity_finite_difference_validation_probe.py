from __future__ import annotations
import json, subprocess, sys, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2248S1198StrictNuBranchGroupPolicySensitivityFiniteDifferenceValidationProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2248_s1198_strict_nu_branch_group_policy_sensitivity_finite_difference_validation_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2248_s1198_strict_nu_branch_group_policy_sensitivity_finite_difference_validation_probe.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(data["schema_version"], "p2248_s1198_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["finite_difference_validation_exported"])
        self.assertTrue(g["validation_pass"])
        self.assertTrue(g["tolerance_strict"])


if __name__ == "__main__":
    unittest.main()
