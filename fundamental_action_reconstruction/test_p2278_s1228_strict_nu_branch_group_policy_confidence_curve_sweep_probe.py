from __future__ import annotations
import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2278S1228StrictNuBranchGroupPolicyConfidenceCurveSweepProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2278_s1228_strict_nu_branch_group_policy_confidence_curve_sweep_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2278_s1228_strict_nu_branch_group_policy_confidence_curve_sweep_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2278_s1228_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["sweep_rows_exported"])
        self.assertTrue(g["all_margins_nonnegative"])
        self.assertTrue(g["all_certified_floors_bounded"])


if __name__ == "__main__":
    unittest.main()
