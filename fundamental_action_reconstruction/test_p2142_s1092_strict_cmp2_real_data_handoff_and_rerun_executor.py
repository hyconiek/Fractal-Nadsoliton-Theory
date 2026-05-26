from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p2142_s1092_strict_cmp2_real_data_handoff_and_rerun_executor.json"


class TestP2142StrictCmp2RealDataHandoffAndRerunExecutor(unittest.TestCase):
    def test_p2142_exports_executor_packet(self) -> None:
        subprocess.run([sys.executable, str(ROOT / "p2142_s1092_strict_cmp2_real_data_handoff_and_rerun_executor.py")], check=True)
        d = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(d["schema_version"], "p2142_s1092_v1")
        self.assertEqual(d["status"], "OPEN_PARTIAL_PROGRESS_WITH_TRACE")
        self.assertIn(d["result_kind"], {
            "PASS_STRICT_CMP2_REAL_DATA_HANDOFF_AND_RERUN_EXECUTOR_WITH_TRACE",
            "OPEN_STRICT_CMP2_REAL_DATA_HANDOFF_AND_RERUN_EXECUTOR_BLOCKED",
        })
        self.assertTrue(d["gatekeeper_checks"]["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
