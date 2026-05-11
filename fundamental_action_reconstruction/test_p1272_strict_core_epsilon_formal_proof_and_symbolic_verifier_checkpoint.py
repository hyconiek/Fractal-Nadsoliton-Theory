#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p1272_strict_core_epsilon_formal_proof_and_symbolic_verifier_checkpoint.py"


class TestP1272ProofAndVerifier(unittest.TestCase):
    def test_outputs_draft_and_partial_verifier(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            p1271 = Path(td) / "p1271.json"
            out = Path(td) / "p1272.json"
            p1271.write_text(json.dumps({"lane": "STRICT_CORE_ONLY"}), encoding="utf-8")
            subprocess.run(["python3", str(SCRIPT), "--p1271", str(p1271), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(payload["formal_proof_text"]["status"], "DRAFT_V1")
            self.assertEqual(payload["symbolic_verifier"]["status"], "PARTIAL_RUN")
            self.assertFalse(payload["strict_kernel_closure_ready"])


if __name__ == "__main__":
    unittest.main()
