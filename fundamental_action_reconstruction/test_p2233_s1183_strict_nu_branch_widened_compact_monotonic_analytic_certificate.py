from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2233S1183StrictNuBranchWidenedCompactMonotonicAnalyticCertificate(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2233_s1183_strict_nu_branch_widened_compact_monotonic_analytic_certificate.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2233_s1183_strict_nu_branch_widened_compact_monotonic_analytic_certificate.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2233_s1183_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["analytic_monotonic_certificate_exported"])
        self.assertTrue(g["analytic_monotonicity_holds"])
        self.assertTrue(g["crossing_persists_on_endpoints"])


if __name__ == "__main__":
    unittest.main()
