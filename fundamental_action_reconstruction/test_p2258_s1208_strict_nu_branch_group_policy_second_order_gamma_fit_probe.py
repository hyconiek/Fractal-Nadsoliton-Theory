from __future__ import annotations
import json, subprocess, sys, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2258S1208StrictNuBranchGroupPolicySecondOrderGammaFitProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2258_s1208_strict_nu_branch_group_policy_second_order_gamma_fit_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2258_s1208_strict_nu_branch_group_policy_second_order_gamma_fit_probe.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(data["schema_version"], "p2258_s1208_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["gamma_fit_exported"])
        self.assertTrue(g["gamma_samples_nonempty"])
        self.assertTrue(g["gamma_spread_nonnegative"])


if __name__ == "__main__":
    unittest.main()
