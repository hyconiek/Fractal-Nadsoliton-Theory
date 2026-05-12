#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p1261_strict_only_uncertainty_decomposition_checkpoint.py"


class TestP1261UncertaintyDecomposition(unittest.TestCase):
    def test_reassesses_all_entries(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            tdp = Path(td)
            p1260 = tdp / "p1260.json"
            out = tdp / "p1261.json"
            p1260.write_text(json.dumps({
                "lane": "STRICT_ONLY",
                "assessments": [
                    {"id": "SP1", "metric": "variance_ratio", "observed": 0.11, "predicted_range": [0.0, 0.15]},
                    {"id": "SP2", "metric": "consistency_pass_rate", "observed": 0.98, "predicted_range": [0.97, 1.0]},
                    {"id": "SP3", "metric": "obstruction_residual_bound", "observed": 0.052, "predicted_range": [0.0, 0.05]},
                ],
            }), encoding="utf-8")

            subprocess.run(["python3", str(SCRIPT), "--p1260", str(p1260), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(payload["lane"], "STRICT_ONLY")
            self.assertEqual(len(payload["reassessed"]), 3)
            self.assertIn("sigma_total", payload["uncertainty_budget"])


if __name__ == "__main__":
    unittest.main()
