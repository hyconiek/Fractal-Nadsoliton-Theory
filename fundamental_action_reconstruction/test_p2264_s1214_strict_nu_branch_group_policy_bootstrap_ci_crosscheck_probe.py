from __future__ import annotations
import json, subprocess, sys, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2264S1214StrictNuBranchGroupPolicyBootstrapCiCrosscheckProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2264_s1214_strict_nu_branch_group_policy_bootstrap_ci_crosscheck_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2264_s1214_strict_nu_branch_group_policy_bootstrap_ci_crosscheck_probe.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(data["schema_version"], "p2264_s1214_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["bootstrap_crosscheck_exported"])
        self.assertTrue(g["bootstrap_ci_ordered"])
        self.assertTrue(g["ci_overlap_computable"])


if __name__ == "__main__":
    unittest.main()
