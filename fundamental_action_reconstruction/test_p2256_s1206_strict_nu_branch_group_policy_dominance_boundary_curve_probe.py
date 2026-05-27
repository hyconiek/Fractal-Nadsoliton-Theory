from __future__ import annotations
import json, subprocess, sys, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2256S1206StrictNuBranchGroupPolicyDominanceBoundaryCurveProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2256_s1206_strict_nu_branch_group_policy_dominance_boundary_curve_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2256_s1206_strict_nu_branch_group_policy_dominance_boundary_curve_probe.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(data["schema_version"], "p2256_s1206_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["boundary_curve_exported"])
        self.assertTrue(g["rows_nonempty"])
        self.assertTrue(g["all_progress_viable_flag_exported"])


if __name__ == "__main__":
    unittest.main()
