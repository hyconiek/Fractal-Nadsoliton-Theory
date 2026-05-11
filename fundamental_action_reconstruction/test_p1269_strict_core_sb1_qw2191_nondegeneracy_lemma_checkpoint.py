#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p1269_strict_core_sb1_qw2191_nondegeneracy_lemma_checkpoint.py"


class TestP1269NonDegeneracy(unittest.TestCase):
    def test_partial_discharge_output(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            p1268 = Path(td) / "p1268.json"
            out = Path(td) / "p1269.json"
            p1268.write_text(json.dumps({"theorem": {"status": "PARTIAL_COMPATIBLE"}}), encoding="utf-8")
            subprocess.run(["python3", str(SCRIPT), "--p1268", str(p1268), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(payload["lemma"]["status"], "PARTIAL_DISCHARGE")
            self.assertFalse(payload["strict_kernel_closure_ready"])


if __name__ == "__main__":
    unittest.main()
