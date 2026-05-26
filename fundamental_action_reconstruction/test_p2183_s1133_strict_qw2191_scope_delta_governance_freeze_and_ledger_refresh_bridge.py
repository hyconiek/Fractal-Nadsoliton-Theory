from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2183S1133StrictQW2191ScopeDeltaGovernanceFreezeAndLedgerRefreshBridge(unittest.TestCase):
    def test_scope_delta_governance_freeze(self) -> None:
        subprocess.run(
            [
                sys.executable,
                str(ROOT / "p2183_s1133_strict_qw2191_scope_delta_governance_freeze_and_ledger_refresh_bridge.py"),
            ],
            check=True,
        )
        d = json.loads(
            (G / "p2183_s1133_strict_qw2191_scope_delta_governance_freeze_and_ledger_refresh_bridge.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(d["schema_version"], "p2183_s1133_v1")
        self.assertTrue(d["gatekeeper_checks"]["scope_delta_governance_freeze_exported"])
        bridge = d["strict_qw2191_scope_delta_governance_freeze_and_ledger_refresh_bridge"]
        self.assertIn("governance_note", bridge)
        self.assertIn("ledger_refresh_bridge_rows", bridge)


if __name__ == "__main__":
    unittest.main()
