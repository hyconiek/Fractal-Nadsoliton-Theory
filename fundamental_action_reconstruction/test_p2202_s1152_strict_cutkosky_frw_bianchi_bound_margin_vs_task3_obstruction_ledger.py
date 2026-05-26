from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2202S1152StrictCutkoskyFrwBianchiBoundMarginVsTask3ObstructionLedger(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2202_s1152_strict_cutkosky_frw_bianchi_bound_margin_vs_task3_obstruction_ledger.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2202_s1152_strict_cutkosky_frw_bianchi_bound_margin_vs_task3_obstruction_ledger.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2202_s1152_v1")
        self.assertTrue(data["gatekeeper_checks"]["obstruction_ledger_exported"])


if __name__ == "__main__":
    unittest.main()
