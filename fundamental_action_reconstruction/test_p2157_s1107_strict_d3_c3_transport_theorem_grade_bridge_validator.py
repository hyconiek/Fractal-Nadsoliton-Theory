from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p2157_s1107_strict_d3_c3_transport_theorem_grade_bridge_validator.json"


class TestP2157StrictD3C3TransportTheoremGradeBridgeValidator(unittest.TestCase):
    def test_export(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2157_s1107_strict_d3_c3_transport_theorem_grade_bridge_validator.py")],
            check=True,
        )
        payload = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(payload["schema_version"], "p2157_s1107_v1")
        self.assertTrue(payload["gatekeeper_checks"]["bridge_validator_exported"])
        self.assertTrue(payload["theorem_grade_bridge_validator"]["theorem_grade_ready"])
        self.assertFalse(payload["gatekeeper_checks"]["full_d3_covariance_transport_proven"])
        self.assertFalse(payload["gatekeeper_checks"]["c3_theorem_proven"])


if __name__ == "__main__":
    unittest.main()
