from __future__ import annotations
import json, subprocess, sys, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2249S1199StrictNuBranchGroupPolicyConservativeStepControllerProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2249_s1199_strict_nu_branch_group_policy_conservative_step_controller_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2249_s1199_strict_nu_branch_group_policy_conservative_step_controller_probe.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(data["schema_version"], "p2249_s1199_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["controller_exported"])
        self.assertTrue(g["base_margin_nonnegative"])
        self.assertTrue(g["guaranteed_margin_nonnegative"])


if __name__ == "__main__":
    unittest.main()
