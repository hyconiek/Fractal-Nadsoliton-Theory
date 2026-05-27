from __future__ import annotations
import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2280S1230StrictNuBranchGroupPolicyMinimalLockCriterionGridProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2280_s1230_strict_nu_branch_group_policy_minimal_lock_criterion_grid_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2280_s1230_strict_nu_branch_group_policy_minimal_lock_criterion_grid_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2280_s1230_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["rows_exported"])
        self.assertTrue(g["cost_proxy_nonnegative"])
        self.assertTrue(g["margins_nonnegative"])
        self.assertTrue(g["certified_floors_bounded"])


if __name__ == "__main__":
    unittest.main()
