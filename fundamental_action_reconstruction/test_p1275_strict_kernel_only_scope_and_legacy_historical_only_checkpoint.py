#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p1275_strict_kernel_only_scope_and_legacy_historical_only_checkpoint.py"


class TestP1275StrictScope(unittest.TestCase):
    def test_strict_only_scope_declared(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            p1274 = Path(td) / "p1274.json"
            out = Path(td) / "p1275.json"
            p1274.write_text(json.dumps({"lane": "STRICT_CORE_ONLY"}), encoding="utf-8")
            subprocess.run(["python3", str(SCRIPT), "--p1274", str(p1274), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(payload["scope_policy"]["active_kernel"], "K_strict_gate")
            self.assertEqual(payload["scope_policy"]["legacy_kernel_operational_use"], "DISALLOWED_IN_CURRENT_STRICT_CORE_PROOF_CHAIN")
            self.assertFalse(payload["strict_kernel_closure_ready"])


if __name__ == "__main__":
    unittest.main()
