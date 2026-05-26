from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2197S1147StrictCutkoskyDriftTrendCadenceAndHardstopPolicy(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2197_s1147_strict_cutkosky_drift_trend_cadence_and_hardstop_policy.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2197_s1147_strict_cutkosky_drift_trend_cadence_and_hardstop_policy.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2197_s1147_v1")
        self.assertTrue(data["gatekeeper_checks"]["cadence_hardstop_policy_exported"])


if __name__ == "__main__":
    unittest.main()
