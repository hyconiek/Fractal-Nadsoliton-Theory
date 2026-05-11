#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p1276_strict_core_local_closure_motion_packet_checkpoint.py"


class TestP1276LocalClosureMotion(unittest.TestCase):
    def test_qw2191_explicitly_not_closed(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            p1275 = Path(td) / "p1275.json"
            p1274 = Path(td) / "p1274.json"
            out = Path(td) / "p1276.json"
            p1275.write_text(json.dumps({"lane": "STRICT_CORE_ONLY"}), encoding="utf-8")
            p1274.write_text(json.dumps({"lane": "STRICT_CORE_ONLY"}), encoding="utf-8")
            subprocess.run(["python3", str(SCRIPT), "--p1275", str(p1275), "--p1274", str(p1274), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(payload["qw2191_closure_status"], "NOT_CLOSED")
            self.assertEqual(payload["global_closure_status"], "OPEN")


if __name__ == "__main__":
    unittest.main()
