from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2214S1164StrictNuBranchClosedFormThresholdConsistencyCertificate(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2214_s1164_strict_nu_branch_closed_form_threshold_consistency_certificate.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2214_s1164_strict_nu_branch_closed_form_threshold_consistency_certificate.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2214_s1164_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["closed_form_consistency_certificate_exported"])
        self.assertTrue(g["closed_form_inside_bisection_bracket"])
        self.assertTrue(g["small_relative_gap_vs_bisection"])


if __name__ == "__main__":
    unittest.main()
