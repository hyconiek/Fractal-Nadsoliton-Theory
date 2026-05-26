from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2200S1150StrictCutkoskyIntegrandBoundAndMonotonicityCertificate(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2200_s1150_strict_cutkosky_integrand_bound_and_monotonicity_certificate.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2200_s1150_strict_cutkosky_integrand_bound_and_monotonicity_certificate.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2200_s1150_v1")
        self.assertTrue(data["gatekeeper_checks"]["bound_certificate_exported"])


if __name__ == "__main__":
    unittest.main()
