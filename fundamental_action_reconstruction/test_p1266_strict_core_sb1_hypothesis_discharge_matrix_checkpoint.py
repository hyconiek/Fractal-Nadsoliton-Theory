#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p1266_strict_core_sb1_hypothesis_discharge_matrix_checkpoint.py"


class TestP1266SB1Matrix(unittest.TestCase):
    def test_emits_open_matrix_for_strict_core(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            p1265 = Path(td) / "p1265.json"
            out = Path(td) / "p1266.json"
            p1265.write_text(json.dumps({"lane": "STRICT_CORE_ONLY"}), encoding="utf-8")
            subprocess.run(["python3", str(SCRIPT), "--p1265", str(p1265), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(payload["lane"], "STRICT_CORE_ONLY")
            self.assertFalse(payload["sb1_discharge_ready"])
            self.assertGreater(payload["open_count"], 0)


if __name__ == "__main__":
    unittest.main()
