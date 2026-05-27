from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2218S1168StrictNuBranchBudgetMapDirectRecomputeMismatchEnvelope(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2218_s1168_strict_nu_branch_budget_map_direct_recompute_mismatch_envelope.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2218_s1168_strict_nu_branch_budget_map_direct_recompute_mismatch_envelope.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2218_s1168_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["mismatch_envelope_exported"])
        self.assertTrue(g["max_abs_mismatch_bounded"])
        self.assertTrue(g["max_rel_mismatch_bounded"])


if __name__ == "__main__":
    unittest.main()
