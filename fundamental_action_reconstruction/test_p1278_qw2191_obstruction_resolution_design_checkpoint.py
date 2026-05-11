#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p1278_qw2191_obstruction_resolution_design_checkpoint.py"


class TestP1278ResolutionDesign(unittest.TestCase):
    def test_emits_three_lemma_repair_plan(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            p1277 = Path(td) / "p1277.json"
            out = Path(td) / "p1278.json"
            p1277.write_text(json.dumps({"qw2191_closure_status": "NOT_CLOSED"}), encoding="utf-8")
            subprocess.run(["python3", str(SCRIPT), "--p1277", str(p1277), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(payload["execution_mode"], "THEOREM_REPAIR_PROGRAM")
            self.assertEqual(len(payload["resolution_design"]["lemmas"]), 3)
            self.assertFalse(payload["strict_kernel_closure_ready"])


if __name__ == "__main__":
    unittest.main()
