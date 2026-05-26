from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p2162_s1112_strict_cmp2_end_to_end_refresh_and_theorem_flag_ledger_snapshot.json"


class TestP2162StrictCmp2EndToEndRefreshAndTheoremFlagLedgerSnapshot(unittest.TestCase):
    def test_export(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2162_s1112_strict_cmp2_end_to_end_refresh_and_theorem_flag_ledger_snapshot.py")],
            check=True,
        )
        payload = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(payload["schema_version"], "p2162_s1112_v1")
        self.assertTrue(payload["gatekeeper_checks"]["ledger_snapshot_exported"])
        self.assertTrue(payload["gatekeeper_checks"]["snapshot_consistent"])


if __name__ == "__main__":
    unittest.main()
