from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p2140_s1090_strict_cmp2_blocked_path_resolution_router.json"


class TestP2140StrictCmp2BlockedPathResolutionRouter(unittest.TestCase):
    def test_p2140_exports_router(self) -> None:
        subprocess.run([sys.executable, str(ROOT / "p2139_s1089_strict_cmp2_nonsynthetic_readiness_freeze_report.py")], check=True)
        subprocess.run([sys.executable, str(ROOT / "p2140_s1090_strict_cmp2_blocked_path_resolution_router.py")], check=True)

        d = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(d["schema_version"], "p2140_s1090_v1")
        self.assertEqual(d["result_kind"], "PASS_STRICT_CMP2_BLOCKED_PATH_RESOLUTION_ROUTER_WITH_TRACE")
        self.assertTrue(d["gatekeeper_checks"]["router_exported"])
        self.assertTrue(d["gatekeeper_checks"]["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
