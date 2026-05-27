from __future__ import annotations
import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2274S1224StrictNuBranchGroupPolicyRobustnessRegionCertificateProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2274_s1224_strict_nu_branch_group_policy_robustness_region_certificate_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2274_s1224_strict_nu_branch_group_policy_robustness_region_certificate_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2274_s1224_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["certified_rows_exported"])
        self.assertTrue(g["all_radius_positive"])
        self.assertTrue(g["all_boxes_within_admissible_domain"])


if __name__ == "__main__":
    unittest.main()
