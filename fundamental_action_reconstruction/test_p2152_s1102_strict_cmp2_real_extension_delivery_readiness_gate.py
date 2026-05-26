from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p2152_s1102_strict_cmp2_real_extension_delivery_readiness_gate.json"


class TestP2152StrictCmp2RealExtensionDeliveryReadinessGate(unittest.TestCase):
    def test_export(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2152_s1102_strict_cmp2_real_extension_delivery_readiness_gate.py")],
            check=True,
        )
        payload = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(payload["schema_version"], "p2152_s1102_v1")
        self.assertTrue(payload["gatekeeper_checks"]["readiness_gate_exported"])


if __name__ == "__main__":
    unittest.main()
