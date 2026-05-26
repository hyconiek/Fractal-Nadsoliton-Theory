from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2193S1143StrictCutkoskyUncertaintyToCiToleranceBridge(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2193_s1143_strict_cutkosky_uncertainty_to_ci_tolerance_bridge.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2193_s1143_strict_cutkosky_uncertainty_to_ci_tolerance_bridge.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2193_s1143_v1")
        self.assertTrue(data["gatekeeper_checks"]["ci_tolerance_bridge_exported"])


if __name__ == "__main__":
    unittest.main()
