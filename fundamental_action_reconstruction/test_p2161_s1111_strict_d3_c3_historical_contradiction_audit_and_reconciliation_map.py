from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p2161_s1111_strict_d3_c3_historical_contradiction_audit_and_reconciliation_map.json"


class TestP2161StrictD3C3HistoricalContradictionAuditAndReconciliationMap(unittest.TestCase):
    def test_export(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2161_s1111_strict_d3_c3_historical_contradiction_audit_and_reconciliation_map.py")],
            check=True,
        )
        payload = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(payload["schema_version"], "p2161_s1111_v1")
        self.assertTrue(payload["gatekeeper_checks"]["reconciliation_map_exported"])
        self.assertTrue(payload["historical_contradiction_audit_and_reconciliation_map"]["contradiction_free_after_reconciliation"])


if __name__ == "__main__":
    unittest.main()
