#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p1277_qw2191_full_closure_theorem_attempt_checkpoint.py"


class TestP1277QW2191ClosureAttempt(unittest.TestCase):
    def test_attempt_reports_obstruction(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            p1276 = Path(td) / "p1276.json"
            out = Path(td) / "p1277.json"
            p1276.write_text(json.dumps({"lane": "STRICT_CORE_ONLY"}), encoding="utf-8")
            subprocess.run(["python3", str(SCRIPT), "--p1276", str(p1276), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(payload["qw2191_closure_status"], "NOT_CLOSED")
            self.assertEqual(payload["theorem_attempt"]["result"], "OBSTRUCTION_REMAINS")
            self.assertEqual(payload["global_closure_status"], "OPEN")


if __name__ == "__main__":
    unittest.main()
