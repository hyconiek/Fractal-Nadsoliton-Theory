from __future__ import annotations
import json, subprocess, sys, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2261S1211StrictNuBranchGroupPolicyStochasticReductionConcentrationBoundProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2261_s1211_strict_nu_branch_group_policy_stochastic_reduction_concentration_bound_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2261_s1211_strict_nu_branch_group_policy_stochastic_reduction_concentration_bound_probe.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(data["schema_version"], "p2261_s1211_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["concentration_bound_exported"])
        self.assertTrue(g["eps_nonnegative"])
        self.assertTrue(g["lower_bound_computable"])


if __name__ == "__main__":
    unittest.main()
