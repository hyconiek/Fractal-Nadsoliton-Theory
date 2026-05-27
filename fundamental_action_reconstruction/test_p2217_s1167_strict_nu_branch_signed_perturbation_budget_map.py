from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2217S1167StrictNuBranchSignedPerturbationBudgetMap(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2217_s1167_strict_nu_branch_signed_perturbation_budget_map.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2217_s1167_strict_nu_branch_signed_perturbation_budget_map.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2217_s1167_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["signed_perturbation_budget_map_exported"])
        self.assertTrue(g["above_budgets_positive_delta"])
        self.assertTrue(g["below_budgets_negative_delta"])


if __name__ == "__main__":
    unittest.main()
