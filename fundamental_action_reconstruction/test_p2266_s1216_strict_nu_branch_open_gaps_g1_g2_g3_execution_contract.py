from __future__ import annotations
import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2266S1216StrictNuBranchOpenGapsG1G2G3ExecutionContract(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2266_s1216_strict_nu_branch_open_gaps_g1_g2_g3_execution_contract.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2266_s1216_strict_nu_branch_open_gaps_g1_g2_g3_execution_contract.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2266_s1216_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["all_three_gaps_contracted"])
        self.assertTrue(g["g1_contracted"])
        self.assertTrue(g["g2_contracted"])
        self.assertTrue(g["g3_contracted"])


if __name__ == "__main__":
    unittest.main()
