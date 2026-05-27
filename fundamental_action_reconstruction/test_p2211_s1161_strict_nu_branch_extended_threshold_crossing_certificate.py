from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2211S1161StrictNuBranchExtendedThresholdCrossingCertificate(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2211_s1161_strict_nu_branch_extended_threshold_crossing_certificate.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2211_s1161_strict_nu_branch_extended_threshold_crossing_certificate.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2211_s1161_v1")
        self.assertTrue(data["gatekeeper_checks"]["extended_crossing_certificate_exported"])
        self.assertTrue(data["gatekeeper_checks"]["crossing_detected_on_extended_grid"])
        self.assertTrue(data["gatekeeper_checks"]["crossing_interval_exported"])


if __name__ == "__main__":
    unittest.main()
