from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2220S1170StrictNuBranchInflatedIntervalCoverageCertificate(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2220_s1170_strict_nu_branch_inflated_interval_coverage_certificate.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2220_s1170_strict_nu_branch_inflated_interval_coverage_certificate.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2220_s1170_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["interval_coverage_certificate_exported"])
        self.assertTrue(g["all_above_targets_covered"])
        self.assertTrue(g["below_coverage_majority"])


if __name__ == "__main__":
    unittest.main()
