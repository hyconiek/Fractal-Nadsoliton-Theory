from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p2138_s1088_strict_cmp2_nonsynthetic_rerun_orchestrator.json"


class TestP2138StrictCmp2NonsyntheticRerunOrchestrator(unittest.TestCase):
    def test_p2138_exports_orchestrator_packet(self) -> None:
        subprocess.run([sys.executable, str(ROOT / "p2138_s1088_strict_cmp2_nonsynthetic_rerun_orchestrator.py")], check=True)

        d = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(d["schema_version"], "p2138_s1088_v1")
        self.assertEqual(d["status"], "OPEN_PARTIAL_PROGRESS_WITH_TRACE")
        self.assertIn(d["result_kind"], {
            "PASS_STRICT_CMP2_NONSYNTHETIC_RERUN_ORCHESTRATOR_WITH_TRACE",
            "OPEN_STRICT_CMP2_NONSYNTHETIC_RERUN_ORCHESTRATOR_BLOCKED",
        })
        self.assertTrue(d["gatekeeper_checks"]["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
