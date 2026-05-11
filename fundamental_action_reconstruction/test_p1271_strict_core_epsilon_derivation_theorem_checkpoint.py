#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p1271_strict_core_epsilon_derivation_theorem_checkpoint.py"


class TestP1271EpsilonDerivation(unittest.TestCase):
    def test_partial_derivation_output(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            p1270 = Path(td) / "p1270.json"
            out = Path(td) / "p1271.json"
            p1270.write_text(json.dumps({"lane": "STRICT_CORE_ONLY"}), encoding="utf-8")
            subprocess.run(["python3", str(SCRIPT), "--p1270", str(p1270), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(payload["theorem"]["status"], "PARTIAL_DISCHARGE")
            self.assertFalse(payload["strict_kernel_closure_ready"])


if __name__ == "__main__":
    unittest.main()
