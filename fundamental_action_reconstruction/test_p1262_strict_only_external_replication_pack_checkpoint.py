#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p1262_strict_only_external_replication_pack_checkpoint.py"
GEN = ROOT / "generated"


class TestP1262ReplicationPack(unittest.TestCase):
    def test_manifest_contains_three_artifacts(self) -> None:
        out = GEN / "_tmp_p1262_test_output.json"
        try:
            subprocess.run(["python3", str(SCRIPT), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(payload["lane"], "STRICT_ONLY")
            self.assertEqual(len(payload["replication_manifest"]), 3)
        finally:
            if out.exists():
                out.unlink()


if __name__ == "__main__":
    unittest.main()
