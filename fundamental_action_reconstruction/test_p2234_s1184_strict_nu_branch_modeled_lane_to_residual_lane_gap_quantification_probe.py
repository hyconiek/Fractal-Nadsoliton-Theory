from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2234S1184StrictNuBranchModeledLaneToResidualLaneGapQuantificationProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2234_s1184_strict_nu_branch_modeled_lane_to_residual_lane_gap_quantification_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2234_s1184_strict_nu_branch_modeled_lane_to_residual_lane_gap_quantification_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2234_s1184_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["gap_quantification_exported"])
        self.assertTrue(g["no_bridge_theorem_claimed"])
        self.assertTrue(g["diagnostic_observables_finite"])


if __name__ == "__main__":
    unittest.main()
