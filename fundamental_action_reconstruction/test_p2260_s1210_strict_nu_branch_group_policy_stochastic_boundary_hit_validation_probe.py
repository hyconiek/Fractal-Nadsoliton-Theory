from __future__ import annotations
import json, subprocess, sys, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2260S1210StrictNuBranchGroupPolicyStochasticBoundaryHitValidationProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2260_s1210_strict_nu_branch_group_policy_stochastic_boundary_hit_validation_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2260_s1210_strict_nu_branch_group_policy_stochastic_boundary_hit_validation_probe.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(data["schema_version"], "p2260_s1210_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["stochastic_validation_exported"])
        self.assertTrue(g["rates_bounded"])
        self.assertTrue(g["reduction_nonnegative"])


if __name__ == "__main__":
    unittest.main()
