from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1284_qw2191_r3_independent_audit_and_replication_checkpoint.py"


class TestP1284QW2191R3IndependentAuditAndReplicationCheckpoint(unittest.TestCase):
    def test_happy_path(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1283 = td_path / "p1283.json"
            out = td_path / "p1284.json"
            p1283.write_text(
                json.dumps(
                    {
                        "next_priority": "R3_INDEPENDENT_AUDIT_AND_REPLICATION",
                        "certificate": {"all_checks_pass": True, "pairwise_checks": [{}, {}, {}]},
                    }
                ),
                encoding="utf-8",
            )
            subprocess.run(["python3", str(SCRIPT), "--p1283", str(p1283), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(payload["independent_audit"]["result"], "PASS")
            self.assertEqual(payload["replication"]["result"], "PASS")
            self.assertEqual(payload["next_priority"], "R4_STRICT_SELECTOR_SOURCE_IDENTIFICATION")

    def test_requires_passed_certificate(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1283 = td_path / "p1283.json"
            out = td_path / "p1284.json"
            p1283.write_text(
                json.dumps(
                    {
                        "next_priority": "R3_INDEPENDENT_AUDIT_AND_REPLICATION",
                        "certificate": {"all_checks_pass": False, "pairwise_checks": [{}, {}]},
                    }
                ),
                encoding="utf-8",
            )
            proc = subprocess.run(["python3", str(SCRIPT), "--p1283", str(p1283), "--out", str(out)], capture_output=True, text=True)
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("all_checks_pass=true", proc.stderr)


if __name__ == "__main__":
    unittest.main()
