from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2208S1158StrictNuBranchSeparationLowerBoundCertificate(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2208_s1158_strict_nu_branch_separation_lower_bound_certificate.py")],
            check=True,
        )
        data = json.loads((G / "p2208_s1158_strict_nu_branch_separation_lower_bound_certificate.json").read_text(encoding="utf-8"))
        self.assertEqual(data["schema_version"], "p2208_s1158_v1")
        self.assertTrue(data["gatekeeper_checks"]["lower_bound_certificate_exported"])
        self.assertTrue(data["gatekeeper_checks"]["linear_bound_holds_on_scan"])


if __name__ == "__main__":
    unittest.main()
