from __future__ import annotations
import json, subprocess, sys, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2268S1218StrictNuBranchGroupPolicyCertificatePassrateLowerBoundProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2268_s1218_strict_nu_branch_group_policy_certificate_passrate_lower_bound_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2268_s1218_strict_nu_branch_group_policy_certificate_passrate_lower_bound_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2268_s1218_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["passrate_lower_bounds_exported"])
        self.assertTrue(g["all_bounds_bounded"])
        self.assertTrue(g["all_targets_bounded"])


if __name__ == "__main__":
    unittest.main()
