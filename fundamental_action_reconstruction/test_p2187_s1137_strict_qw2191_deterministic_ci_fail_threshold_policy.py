from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2187S1137StrictQW2191DeterministicCiFailThresholdPolicy(unittest.TestCase):
    def test_policy_export_and_decision_fields(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2187_s1137_strict_qw2191_deterministic_ci_fail_threshold_policy.py")],
            check=True,
        )
        data = json.loads(
            (
                G / "p2187_s1137_strict_qw2191_deterministic_ci_fail_threshold_policy.json"
            ).read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2187_s1137_v1")
        self.assertTrue(data["gatekeeper_checks"]["deterministic_policy_exported"])
        policy = data["strict_qw2191_deterministic_ci_fail_threshold_policy"]
        self.assertIn(policy["final_ci_decision"], {"PASS", "WARN", "FAIL"})
        self.assertIn("policy_table", policy)


if __name__ == "__main__":
    unittest.main()
