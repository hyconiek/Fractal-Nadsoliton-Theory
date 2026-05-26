from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2181S1131StrictQW2191ReplayCertificateFreeze(unittest.TestCase):
    def test_replay_certificate(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2181_s1131_strict_qw2191_replay_certificate_freeze.py")],
            check=True,
        )
        d = json.loads((G / "p2181_s1131_strict_qw2191_replay_certificate_freeze.json").read_text(encoding="utf-8"))
        self.assertEqual(d["schema_version"], "p2181_s1131_v1")
        self.assertTrue(d["gatekeeper_checks"]["replay_certificate_exported"])
        cert = d["strict_qw2191_replay_certificate_freeze"]["replay_certificate"]
        self.assertIn("certificate_id", cert)


if __name__ == "__main__":
    unittest.main()
