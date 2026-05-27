from __future__ import annotations
import json, subprocess, sys, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2262S1212StrictNuBranchGroupPolicyMultiseedReductionStabilityProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2262_s1212_strict_nu_branch_group_policy_multiseed_reduction_stability_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2262_s1212_strict_nu_branch_group_policy_multiseed_reduction_stability_probe.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(data["schema_version"], "p2262_s1212_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["multiseed_scan_exported"])
        self.assertTrue(g["all_seed_reductions_nonnegative"])
        self.assertTrue(g["reduction_interval_ordered"])


if __name__ == "__main__":
    unittest.main()
