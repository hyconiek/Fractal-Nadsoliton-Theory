#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p1273_strict_core_d2_symbolic_completion_and_proof_promotion_checkpoint.py"


class TestP1273D2Completion(unittest.TestCase):
    def test_promotes_proof_and_keeps_global_closure_blocked(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            p1272 = Path(td) / "p1272.json"
            out = Path(td) / "p1273.json"
            p1272.write_text(json.dumps({"lane": "STRICT_CORE_ONLY"}), encoding="utf-8")
            subprocess.run(["python3", str(SCRIPT), "--p1272", str(p1272), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(payload["proof_promotion"]["to"], "REVIEWED_FORMAL_PROOF")
            self.assertEqual(payload["d2_symbolic_verification"]["status"], "COMPLETE")
            self.assertFalse(payload["strict_kernel_closure_ready"])


if __name__ == "__main__":
    unittest.main()
