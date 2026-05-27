from __future__ import annotations
import json, subprocess, sys, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2254S1204StrictNuBranchGroupPolicyAdaptiveRhoDominanceConditionProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2254_s1204_strict_nu_branch_group_policy_adaptive_rho_dominance_condition_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2254_s1204_strict_nu_branch_group_policy_adaptive_rho_dominance_condition_probe.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(data["schema_version"], "p2254_s1204_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["dominance_condition_exported"])
        self.assertTrue(g["progress_gap_computable"])
        self.assertTrue(g["margin_floor_gap_computable"])


if __name__ == "__main__":
    unittest.main()
