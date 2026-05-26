from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2195S1145StrictCutkoskyMultiBackgroundReleaseGateOverridePolicy(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2195_s1145_strict_cutkosky_multi_background_release_gate_override_policy.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2195_s1145_strict_cutkosky_multi_background_release_gate_override_policy.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2195_s1145_v1")
        self.assertTrue(data["gatekeeper_checks"]["override_policy_exported"])


if __name__ == "__main__":
    unittest.main()
