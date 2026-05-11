#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p1265_strict_core_symmetry_breaking_source_theorem_scaffold.py"


class TestP1265SB1Scaffold(unittest.TestCase):
    def test_emits_strict_core_scaffold(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            out = Path(td) / "p1265.json"
            subprocess.run(["python3", str(SCRIPT), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(payload["lane"], "STRICT_CORE_ONLY")
            self.assertEqual(payload["theorem"]["id"], "SB1")
            self.assertEqual(payload["qw2191_interface"]["status"], "OPEN")


if __name__ == "__main__":
    unittest.main()
