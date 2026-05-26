from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p2143_s1093_strict_cmp2_real_ci_stability_readiness_memo.json"


class TestP2143StrictCmp2RealCiStabilityReadinessMemo(unittest.TestCase):
    def test_p2143_exports_readiness_memo(self) -> None:
        subprocess.run([sys.executable, str(ROOT / "p2143_s1093_strict_cmp2_real_ci_stability_readiness_memo.py")], check=True)
        d = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(d["schema_version"], "p2143_s1093_v1")
        self.assertIn(d["result_kind"], {
            "PASS_STRICT_CMP2_REAL_CI_STABILITY_READINESS_MEMO_WITH_TRACE",
            "OPEN_STRICT_CMP2_REAL_CI_STABILITY_READINESS_MEMO_BLOCKED",
        })
        self.assertTrue(d["gatekeeper_checks"]["readiness_memo_exported"])


if __name__ == "__main__":
    unittest.main()
