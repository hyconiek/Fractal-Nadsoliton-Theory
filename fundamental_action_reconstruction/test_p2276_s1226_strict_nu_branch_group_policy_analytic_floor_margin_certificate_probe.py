from __future__ import annotations
import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2276S1226StrictNuBranchGroupPolicyAnalyticFloorMarginCertificateProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2276_s1226_strict_nu_branch_group_policy_analytic_floor_margin_certificate_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2276_s1226_strict_nu_branch_group_policy_analytic_floor_margin_certificate_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2276_s1226_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["certificate_rows_exported"])
        self.assertTrue(g["all_targets_bounded"])
        self.assertTrue(g["all_certified_floors_bounded"])
        self.assertTrue(g["confidence_margin_nonnegative"])


if __name__ == "__main__":
    unittest.main()
