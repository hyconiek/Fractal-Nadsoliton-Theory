from __future__ import annotations
import json, subprocess, sys, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2247S1197StrictNuBranchGroupPolicyBoundarySensitivityProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2247_s1197_strict_nu_branch_group_policy_boundary_sensitivity_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2247_s1197_strict_nu_branch_group_policy_boundary_sensitivity_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2247_s1197_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["boundary_sensitivity_exported"])
        self.assertTrue(g["derivative_wrt_group_count_negative"])
        self.assertTrue(g["derivative_wrt_load_negative"])


if __name__ == "__main__":
    unittest.main()
