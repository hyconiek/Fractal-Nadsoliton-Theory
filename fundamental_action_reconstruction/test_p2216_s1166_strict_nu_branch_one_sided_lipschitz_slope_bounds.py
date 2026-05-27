from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2216S1166StrictNuBranchOneSidedLipschitzSlopeBounds(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2216_s1166_strict_nu_branch_one_sided_lipschitz_slope_bounds.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2216_s1166_strict_nu_branch_one_sided_lipschitz_slope_bounds.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2216_s1166_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["one_sided_slope_bounds_exported"])
        self.assertTrue(g["left_local_slope_positive"])
        self.assertTrue(g["right_local_slope_positive"])


if __name__ == "__main__":
    unittest.main()
