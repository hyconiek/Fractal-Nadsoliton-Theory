#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p1258_strict_only_operational_lane_commitment_checkpoint.py"


class TestP1258StrictOnlyLane(unittest.TestCase):
    def test_outputs_strict_only_with_hard_constraints(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            tdp = Path(td)
            p1256 = tdp / "p1256.json"
            out = tdp / "p1258.json"
            p1256.write_text(json.dumps({
                "closure_policy": "STRICT_CLOSURE_FORBIDDEN_UNTIL_B1_OR_NB1_DISCHARGED",
            }), encoding="utf-8")

            subprocess.run(["python3", str(SCRIPT), "--p1256", str(p1256), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))

            self.assertEqual(payload["strategy_decision"]["operational_lane"], "STRICT_ONLY")
            self.assertIn("NO_STRICT_CORE_CLOSURE_UNTIL_B1_OR_NB1", payload["hard_constraints"])


if __name__ == "__main__":
    unittest.main()
