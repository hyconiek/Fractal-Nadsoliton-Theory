from __future__ import annotations
import json, subprocess, sys, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2240S1190StrictNuBranchGroupPolicySignStabilityConditionProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2240_s1190_strict_nu_branch_group_policy_sign_stability_condition_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2240_s1190_strict_nu_branch_group_policy_sign_stability_condition_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2240_s1190_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["group_sign_stability_condition_exported"])
        self.assertTrue(g["all_groups_satisfy_sufficient_condition"])
        self.assertTrue(g["signed_margin_nonnegative"])


if __name__ == "__main__":
    unittest.main()
