from __future__ import annotations
import json, subprocess, sys, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2259S1209StrictNuBranchGroupPolicyGammaInjectedBoundaryHitRateProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2259_s1209_strict_nu_branch_group_policy_gamma_injected_boundary_hit_rate_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2259_s1209_strict_nu_branch_group_policy_gamma_injected_boundary_hit_rate_probe.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(data["schema_version"], "p2259_s1209_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["gamma_injected_hit_rate_exported"])
        self.assertTrue(g["hit_rates_bounded"])
        self.assertTrue(g["hit_rate_reduction_computable"])


if __name__ == "__main__":
    unittest.main()
