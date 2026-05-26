from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2194S1144StrictCutkoskyMultiBackgroundToleranceStabilityAudit(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2194_s1144_strict_cutkosky_multi_background_tolerance_stability_audit.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2194_s1144_strict_cutkosky_multi_background_tolerance_stability_audit.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2194_s1144_v1")
        self.assertTrue(data["gatekeeper_checks"]["multi_background_audit_exported"])


if __name__ == "__main__":
    unittest.main()
