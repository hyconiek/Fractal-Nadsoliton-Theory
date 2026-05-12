#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p1255_bridge_or_nonbridge_decision_gate_checkpoint.py"


class TestP1255BridgeDecisionGate(unittest.TestCase):
    def test_gate_passes_when_bridge_unresolved_and_closure_open(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            in_path = td_path / "p1254.json"
            out_path = td_path / "p1255.json"

            in_payload = {
                "theory_closure_status": "OPEN",
                "strict_closure_claim_allowed": False,
                "open_obligation_count_after_l2_attempt": 1,
            }
            in_path.write_text(json.dumps(in_payload), encoding="utf-8")

            subprocess.run([
                "python3", str(SCRIPT), "--p1254", str(in_path), "--out", str(out_path)
            ], check=True)

            out = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertTrue(out["gate_pass"])
            self.assertEqual(out["decision"], "PROCEED_TO_BRIDGE_OR_NONBRIDGE_FORMALIZATION_ONLY")
            self.assertEqual(out["bridge_status"], "UNRESOLVED")


if __name__ == "__main__":
    unittest.main()
