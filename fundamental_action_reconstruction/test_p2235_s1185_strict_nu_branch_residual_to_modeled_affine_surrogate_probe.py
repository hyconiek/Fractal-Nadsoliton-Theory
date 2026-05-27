from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2235S1185StrictNuBranchResidualToModeledAffineSurrogateProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2235_s1185_strict_nu_branch_residual_to_modeled_affine_surrogate_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2235_s1185_strict_nu_branch_residual_to_modeled_affine_surrogate_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2235_s1185_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["affine_surrogate_exported"])
        self.assertTrue(g["positive_slope_map"])
        self.assertTrue(g["all_surrogates_inside_modeled_range"])


if __name__ == "__main__":
    unittest.main()
