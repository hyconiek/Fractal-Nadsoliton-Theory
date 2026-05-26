from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2179S1129StrictQW2191CompactLedgerFreezeAndAuditHandoff(unittest.TestCase):
    def test_compact_ledger_freeze(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2179_s1129_strict_qw2191_compact_ledger_freeze_and_audit_handoff.py")],
            check=True,
        )
        d = json.loads(
            (G / "p2179_s1129_strict_qw2191_compact_ledger_freeze_and_audit_handoff.json").read_text(encoding="utf-8")
        )
        self.assertEqual(d["schema_version"], "p2179_s1129_v1")
        self.assertTrue(d["gatekeeper_checks"]["compact_ledger_freeze_exported"])
        handoff = d["strict_qw2191_compact_ledger_freeze_and_audit_handoff"]["audit_handoff"]
        self.assertIn("handoff_contract_id", handoff)


if __name__ == "__main__":
    unittest.main()
