#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p1259_strict_only_prediction_risk_ledger_checkpoint.py"


class TestP1259StrictPredictionLedger(unittest.TestCase):
    def test_emits_falsifiable_predictions(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            tdp = Path(td)
            p1258 = tdp / "p1258.json"
            out = tdp / "p1259.json"
            p1258.write_text(json.dumps({"strategy_decision": {"operational_lane": "STRICT_ONLY"}}), encoding="utf-8")

            subprocess.run(["python3", str(SCRIPT), "--p1258", str(p1258), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))

            self.assertEqual(payload["lane"], "STRICT_ONLY")
            self.assertEqual(len(payload["predictions"]), 3)
            self.assertTrue(all(p["status"] == "OPEN_TEST" for p in payload["predictions"]))


if __name__ == "__main__":
    unittest.main()
