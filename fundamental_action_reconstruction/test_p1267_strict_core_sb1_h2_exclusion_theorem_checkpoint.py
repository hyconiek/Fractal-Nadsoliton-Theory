#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p1267_strict_core_sb1_h2_exclusion_theorem_checkpoint.py"


class TestP1267H2Exclusion(unittest.TestCase):
    def test_h2_partial_discharge_output(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            p1266 = Path(td) / "p1266.json"
            out = Path(td) / "p1267.json"
            p1266.write_text(json.dumps({
                "sb1_hypothesis_matrix": [{"hypothesis": "H2", "status": "OPEN"}],
            }), encoding="utf-8")
            subprocess.run(["python3", str(SCRIPT), "--p1266", str(p1266), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(payload["lane"], "STRICT_CORE_ONLY")
            self.assertEqual(payload["h2_status_update"], "PARTIAL")
            self.assertFalse(payload["strict_kernel_closure_ready"])


if __name__ == "__main__":
    unittest.main()
