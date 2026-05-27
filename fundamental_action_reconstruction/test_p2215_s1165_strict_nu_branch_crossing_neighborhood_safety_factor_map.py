from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2215S1165StrictNuBranchCrossingNeighborhoodSafetyFactorMap(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2215_s1165_strict_nu_branch_crossing_neighborhood_safety_factor_map.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2215_s1165_strict_nu_branch_crossing_neighborhood_safety_factor_map.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2215_s1165_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["safety_factor_map_exported"])
        self.assertTrue(g["below_dm_star_not_certified"])
        self.assertTrue(g["above_dm_star_certified"])


if __name__ == "__main__":
    unittest.main()
