from __future__ import annotations
import json, subprocess, sys, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2244S1194StrictNuBranchGroupPolicyViolationTailBoundProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2244_s1194_strict_nu_branch_group_policy_violation_tail_bound_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2244_s1194_strict_nu_branch_group_policy_violation_tail_bound_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2244_s1194_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["tail_bound_exported"])
        self.assertTrue(g["deterministic_bound_computable"])
        self.assertTrue(g["statistical_bound_computable"])


if __name__ == "__main__":
    unittest.main()
