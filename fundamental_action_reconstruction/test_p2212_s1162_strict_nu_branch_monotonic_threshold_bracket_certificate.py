from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2212S1162StrictNuBranchMonotonicThresholdBracketCertificate(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2212_s1162_strict_nu_branch_monotonic_threshold_bracket_certificate.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2212_s1162_strict_nu_branch_monotonic_threshold_bracket_certificate.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2212_s1162_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["monotonic_bracket_certificate_exported"])
        self.assertTrue(g["monotonicity_verified_on_bracket"])
        self.assertTrue(g["lo_uncertified_endpoint_verified"])
        self.assertTrue(g["hi_certified_endpoint_verified"])


if __name__ == "__main__":
    unittest.main()
