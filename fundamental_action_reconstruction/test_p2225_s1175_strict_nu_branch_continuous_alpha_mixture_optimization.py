from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2225S1175StrictNuBranchContinuousAlphaMixtureOptimization(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2225_s1175_strict_nu_branch_continuous_alpha_mixture_optimization.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2225_s1175_strict_nu_branch_continuous_alpha_mixture_optimization.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2225_s1175_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["continuous_alpha_optimization_exported"])
        self.assertTrue(g["alpha_star_in_unit_interval"])
        self.assertTrue(g["worst_deficit_zero_at_alpha_star"])


if __name__ == "__main__":
    unittest.main()
