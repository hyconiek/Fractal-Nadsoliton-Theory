from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2180S1130StrictQW2191IndependentAuditReplay(unittest.TestCase):
    def test_independent_audit_replay(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2180_s1130_strict_qw2191_independent_audit_replay_against_handoff_contract.py")],
            check=True,
        )
        d = json.loads(
            (G / "p2180_s1130_strict_qw2191_independent_audit_replay_against_handoff_contract.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(d["schema_version"], "p2180_s1130_v1")
        self.assertTrue(d["gatekeeper_checks"]["independent_audit_replay_exported"])
        self.assertIn(
            "checks",
            d["strict_qw2191_independent_audit_replay_against_handoff_contract"],
        )


if __name__ == "__main__":
    unittest.main()
