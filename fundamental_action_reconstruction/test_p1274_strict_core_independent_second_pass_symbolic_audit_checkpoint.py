#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p1274_strict_core_independent_second_pass_symbolic_audit_checkpoint.py"


class TestP1274SecondPassAudit(unittest.TestCase):
    def test_audit_pass_and_global_gate_still_closed(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            p1273 = Path(td) / "p1273.json"
            out = Path(td) / "p1274.json"
            p1273.write_text(json.dumps({"lane": "STRICT_CORE_ONLY"}), encoding="utf-8")
            subprocess.run(["python3", str(SCRIPT), "--p1273", str(p1273), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(payload["audit"]["status"], "PASS")
            self.assertTrue(payload["audit"]["consistency_with_pass1"])
            self.assertFalse(payload["strict_kernel_closure_ready"])


if __name__ == "__main__":
    unittest.main()
