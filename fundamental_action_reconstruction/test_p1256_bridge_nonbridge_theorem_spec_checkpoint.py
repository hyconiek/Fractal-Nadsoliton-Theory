#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p1256_bridge_nonbridge_theorem_spec_checkpoint.py"


class TestP1256BridgeNonbridgeSpec(unittest.TestCase):
    def test_emits_open_obligations_and_closure_policy(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1255 = td_path / "p1255.json"
            out = td_path / "p1256.json"
            p1255.write_text(json.dumps({
                "decision": "PROCEED_TO_BRIDGE_OR_NONBRIDGE_FORMALIZATION_ONLY",
                "gate_pass": True,
            }), encoding="utf-8")

            subprocess.run(["python3", str(SCRIPT), "--p1255", str(p1255), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))

            self.assertEqual(payload["packet"], "P1256")
            self.assertEqual(payload["closure_policy"], "STRICT_CLOSURE_FORBIDDEN_UNTIL_B1_OR_NB1_DISCHARGED")
            self.assertTrue(all(o["status"] == "OPEN" for o in payload["proof_obligations"]))


if __name__ == "__main__":
    unittest.main()
