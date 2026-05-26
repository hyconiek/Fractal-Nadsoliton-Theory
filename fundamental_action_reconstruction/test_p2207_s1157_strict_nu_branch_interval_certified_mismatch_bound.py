from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2207S1157StrictNuBranchIntervalCertifiedMismatchBound(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2207_s1157_strict_nu_branch_interval_certified_mismatch_bound.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2207_s1157_strict_nu_branch_interval_certified_mismatch_bound.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2207_s1157_v1")
        self.assertTrue(data["gatekeeper_checks"]["interval_mismatch_certificate_exported"])
        self.assertTrue(data["gatekeeper_checks"]["nontrivial_branch_separation_positive"])


if __name__ == "__main__":
    unittest.main()
