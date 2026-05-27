from __future__ import annotations
import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2275S1225StrictNuBranchGroupPolicyCertifiedBoxCornerReplayPassrateFloorProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2275_s1225_strict_nu_branch_group_policy_certified_box_corner_replay_passrate_floor_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2275_s1225_strict_nu_branch_group_policy_certified_box_corner_replay_passrate_floor_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2275_s1225_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["replay_rows_exported"])
        self.assertTrue(g["all_corner_stats_nonempty"])
        self.assertTrue(g["all_target_floors_bounded"])
        self.assertTrue(g["all_empirical_floors_bounded"])


if __name__ == "__main__":
    unittest.main()
