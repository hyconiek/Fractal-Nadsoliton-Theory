#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p1268_strict_core_sb1_qw2191_compatibility_theorem_checkpoint.py"


class TestP1268Compatibility(unittest.TestCase):
    def test_partial_compatible_output(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            p1266 = Path(td) / "p1266.json"
            p1267 = Path(td) / "p1267.json"
            out = Path(td) / "p1268.json"
            p1266.write_text(json.dumps({"lane": "STRICT_CORE_ONLY"}), encoding="utf-8")
            p1267.write_text(json.dumps({"lane": "STRICT_CORE_ONLY"}), encoding="utf-8")
            subprocess.run(["python3", str(SCRIPT), "--p1266", str(p1266), "--p1267", str(p1267), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(payload["theorem"]["status"], "PARTIAL_COMPATIBLE")
            self.assertFalse(payload["strict_kernel_closure_ready"])


if __name__ == "__main__":
    unittest.main()
