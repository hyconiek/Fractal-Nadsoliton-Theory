from __future__ import annotations
import json, subprocess, sys, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2269S1219StrictNuBranchGroupPolicyPassrateBoundMonotonicityLipschitzProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2269_s1219_strict_nu_branch_group_policy_passrate_bound_monotonicity_lipschitz_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2269_s1219_strict_nu_branch_group_policy_passrate_bound_monotonicity_lipschitz_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2269_s1219_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["monotonicity_lipschitz_exported"])
        self.assertTrue(g["all_rows_monotone_in_rho"])
        self.assertTrue(g["all_rows_monotone_in_kappa"])
        self.assertTrue(g["lipschitz_l1_nonnegative"])


if __name__ == "__main__":
    unittest.main()
