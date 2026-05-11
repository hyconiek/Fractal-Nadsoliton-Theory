#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p1260_strict_only_benchmark_execution_checkpoint.py"


class TestP1260BenchmarkExecution(unittest.TestCase):
    def test_emits_assessments_and_counts(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            tdp = Path(td)
            p1259 = tdp / "p1259.json"
            out = tdp / "p1260.json"
            p1259.write_text(json.dumps({
                "lane": "STRICT_ONLY",
                "predictions": [
                    {"id": "SP1", "metric": "variance_ratio", "predicted_range": [0.0, 0.15], "falsification_rule": "x"},
                    {"id": "SP2", "metric": "consistency_pass_rate", "predicted_range": [0.97, 1.0], "falsification_rule": "x"},
                    {"id": "SP3", "metric": "obstruction_residual_bound", "predicted_range": [0.0, 0.05], "falsification_rule": "x"},
                ],
            }), encoding="utf-8")

            subprocess.run(["python3", str(SCRIPT), "--p1259", str(p1259), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))

            self.assertEqual(payload["lane"], "STRICT_ONLY")
            self.assertEqual(len(payload["assessments"]), 3)
            self.assertIn("INCONCLUSIVE", payload["status_counts"])


if __name__ == "__main__":
    unittest.main()
