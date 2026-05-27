from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2232S1182StrictNuBranchWidenedCompactIntervalMonotoneExtrapolationProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2232_s1182_strict_nu_branch_widened_compact_interval_monotone_extrapolation_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2232_s1182_strict_nu_branch_widened_compact_interval_monotone_extrapolation_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2232_s1182_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["widened_compact_probe_exported"])
        self.assertTrue(g["monotone_on_widened_interval"])
        self.assertTrue(g["crossing_detected_on_widened_interval"])


if __name__ == "__main__":
    unittest.main()
