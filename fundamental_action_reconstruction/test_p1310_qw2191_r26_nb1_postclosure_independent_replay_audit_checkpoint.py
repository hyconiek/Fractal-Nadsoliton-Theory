from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1310_qw2191_r26_nb1_postclosure_independent_replay_audit_checkpoint.py"


class TestP1310QW2191R26NB1PostclosureIndependentReplayAuditCheckpoint(unittest.TestCase):
    def test_happy_path(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1309 = td_path / "p1309.json"
            out = td_path / "p1310.json"
            p1309.write_text(
                json.dumps(
                    {
                        "next_priority": "R26_NB1_POSTCLOSURE_INDEPENDENT_REPLAY_AUDIT",
                        "r25_closure_statement": {"status": "FORMAL_CLOSURE_DECLARED"},
                    }
                ),
                encoding="utf-8",
            )
            subprocess.run(["python3", str(SCRIPT), "--p1309", str(p1309), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(payload["r26_independent_replay_audit"]["status"], "AUDIT_COMPLETE")

    def test_requires_closure_declared(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1309 = td_path / "p1309.json"
            out = td_path / "p1310.json"
            p1309.write_text(
                json.dumps(
                    {
                        "next_priority": "R26_NB1_POSTCLOSURE_INDEPENDENT_REPLAY_AUDIT",
                        "r25_closure_statement": {"status": "BLOCKED"},
                    }
                ),
                encoding="utf-8",
            )
            proc = subprocess.run(["python3", str(SCRIPT), "--p1309", str(p1309), "--out", str(out)], capture_output=True, text=True)
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("FORMAL_CLOSURE_DECLARED", proc.stderr)


if __name__ == "__main__":
    unittest.main()
