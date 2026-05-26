from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2189S1139StrictQW2191PolicyReleaseContractReplayInvarianceAudit(unittest.TestCase):
    def test_packet_and_invariance_checks(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2189_s1139_strict_qw2191_policy_release_contract_replay_invariance_audit.py")],
            check=True,
        )
        data = json.loads(
            (
                G / "p2189_s1139_strict_qw2191_policy_release_contract_replay_invariance_audit.json"
            ).read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2189_s1139_v1")
        self.assertTrue(data["gatekeeper_checks"]["replay_invariance_audit_exported"])
        self.assertIn("strict_qw2191_policy_release_contract_replay_invariance_audit", data)


if __name__ == "__main__":
    unittest.main()
