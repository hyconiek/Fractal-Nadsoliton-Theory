from __future__ import annotations
import json, subprocess, sys, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2257S1207StrictNuBranchGroupPolicySecondOrderBoundaryShiftProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2257_s1207_strict_nu_branch_group_policy_second_order_boundary_shift_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2257_s1207_strict_nu_branch_group_policy_second_order_boundary_shift_probe.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(data["schema_version"], "p2257_s1207_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["second_order_boundary_exported"])
        self.assertTrue(g["rows_nonempty"])
        self.assertTrue(g["max_shift_nonnegative"])


if __name__ == "__main__":
    unittest.main()
