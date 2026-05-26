from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2196S1146StrictCutkoskyDriftFrequencyTrendReplay(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2196_s1146_strict_cutkosky_drift_frequency_trend_replay.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2196_s1146_strict_cutkosky_drift_frequency_trend_replay.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2196_s1146_v1")
        self.assertTrue(data["gatekeeper_checks"]["drift_frequency_trend_exported"])


if __name__ == "__main__":
    unittest.main()
