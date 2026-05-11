#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p1270_strict_core_qw2191_full_obstruction_bound_proof_checkpoint.py"


class TestP1270FullObstructionBound(unittest.TestCase):
    def test_partial_bound_proof_output(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            p1269 = Path(td) / "p1269.json"
            out = Path(td) / "p1270.json"
            p1269.write_text(json.dumps({"lane": "STRICT_CORE_ONLY"}), encoding="utf-8")
            subprocess.run(["python3", str(SCRIPT), "--p1269", str(p1269), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(payload["bound_proof"]["status"], "PARTIAL_DISCHARGE")
            self.assertFalse(payload["strict_kernel_closure_ready"])


if __name__ == "__main__":
    unittest.main()
