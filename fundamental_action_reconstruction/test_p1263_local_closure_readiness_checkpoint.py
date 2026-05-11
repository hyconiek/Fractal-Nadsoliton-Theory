#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p1263_local_closure_readiness_checkpoint.py"


class TestP1263LocalClosureReadiness(unittest.TestCase):
    def test_emits_local_ready_and_global_open(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            out = Path(td) / "p1263.json"
            subprocess.run(["python3", str(SCRIPT), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(payload["local_readiness_status"], "READY_FOR_LOCAL_HANDOFF")
            self.assertEqual(payload["global_theory_closure_status"], "OPEN")


if __name__ == "__main__":
    unittest.main()
