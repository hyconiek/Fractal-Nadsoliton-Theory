from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2213S1163StrictNuBranchBisectionThresholdToleranceCertificate(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2213_s1163_strict_nu_branch_bisection_threshold_tolerance_certificate.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2213_s1163_strict_nu_branch_bisection_threshold_tolerance_certificate.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2213_s1163_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["bisection_certificate_exported"])
        self.assertTrue(g["input_bracket_straddles_crossing"])
        self.assertTrue(g["final_width_within_tolerance"])


if __name__ == "__main__":
    unittest.main()
