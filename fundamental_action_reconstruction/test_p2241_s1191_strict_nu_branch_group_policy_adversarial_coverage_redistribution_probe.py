from __future__ import annotations
import json, subprocess, sys, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2241S1191StrictNuBranchGroupPolicyAdversarialCoverageRedistributionProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2241_s1191_strict_nu_branch_group_policy_adversarial_coverage_redistribution_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2241_s1191_strict_nu_branch_group_policy_adversarial_coverage_redistribution_probe.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(data["schema_version"], "p2241_s1191_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["adversarial_redistribution_exported"])
        self.assertTrue(g["baseline_condition_holds"])
        self.assertIn("adversarial_fragility_detected", g)
        self.assertIn("adversarial_condition_holds", g)


if __name__ == "__main__":
    unittest.main()
