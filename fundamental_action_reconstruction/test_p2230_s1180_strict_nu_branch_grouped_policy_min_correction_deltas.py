from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2230S1180StrictNuBranchGroupedPolicyMinCorrectionDeltas(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2230_s1180_strict_nu_branch_grouped_policy_min_correction_deltas.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2230_s1180_strict_nu_branch_grouped_policy_min_correction_deltas.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2230_s1180_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["correction_delta_certificate_exported"])
        self.assertTrue(g["deltas_nonnegative"])
        self.assertTrue(g["corrected_policy_meets_required"])


if __name__ == "__main__":
    unittest.main()
