from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2210S1160StrictNuBranchThresholdFlipIntervalCertificate(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2210_s1160_strict_nu_branch_threshold_flip_interval_certificate.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2210_s1160_strict_nu_branch_threshold_flip_interval_certificate.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2210_s1160_v1")
        self.assertTrue(data["gatekeeper_checks"]["flip_interval_certificate_exported"])
        self.assertTrue(data["gatekeeper_checks"]["coarse_rows_present"])
        self.assertTrue(data["gatekeeper_checks"]["no_flip_detected_in_coarse_window"])


if __name__ == "__main__":
    unittest.main()
