from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2229S1179StrictNuBranchGroupedPolicyMissSetCertificate(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2229_s1179_strict_nu_branch_grouped_policy_miss_set_certificate.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2229_s1179_strict_nu_branch_grouped_policy_miss_set_certificate.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2229_s1179_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["miss_set_certificate_exported"])
        self.assertTrue(g["all_rows_covered"])
        self.assertTrue(g["miss_set_empty"])


if __name__ == "__main__":
    unittest.main()
