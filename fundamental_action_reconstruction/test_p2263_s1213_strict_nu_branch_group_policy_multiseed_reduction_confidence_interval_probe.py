from __future__ import annotations
import json, subprocess, sys, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2263S1213StrictNuBranchGroupPolicyMultiseedReductionConfidenceIntervalProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2263_s1213_strict_nu_branch_group_policy_multiseed_reduction_confidence_interval_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2263_s1213_strict_nu_branch_group_policy_multiseed_reduction_confidence_interval_probe.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(data["schema_version"], "p2263_s1213_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["confidence_interval_exported"])
        self.assertTrue(g["ci_bounds_ordered"])
        self.assertTrue(g["std_nonnegative"])


if __name__ == "__main__":
    unittest.main()
