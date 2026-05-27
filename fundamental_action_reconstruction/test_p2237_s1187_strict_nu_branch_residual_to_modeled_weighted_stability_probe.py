from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2237S1187StrictNuBranchResidualToModeledWeightedStabilityProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2237_s1187_strict_nu_branch_residual_to_modeled_weighted_stability_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2237_s1187_strict_nu_branch_residual_to_modeled_weighted_stability_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2237_s1187_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["weighted_stability_probe_exported"])
        self.assertTrue(g["weighted_slope_positive"])
        self.assertTrue(g["train_split_slope_positive"])
        self.assertTrue(g["test_points_inside_train_band"])


if __name__ == "__main__":
    unittest.main()
