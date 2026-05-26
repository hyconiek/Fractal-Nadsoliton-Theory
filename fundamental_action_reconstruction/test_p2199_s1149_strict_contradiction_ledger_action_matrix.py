from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2199S1149StrictContradictionLedgerActionMatrix(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2199_s1149_strict_contradiction_ledger_action_matrix.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2199_s1149_strict_contradiction_ledger_action_matrix.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2199_s1149_v1")
        self.assertTrue(data["gatekeeper_checks"]["action_matrix_exported"])


if __name__ == "__main__":
    unittest.main()
