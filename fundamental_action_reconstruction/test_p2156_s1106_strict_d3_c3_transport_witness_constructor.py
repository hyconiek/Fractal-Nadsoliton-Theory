from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p2156_s1106_strict_d3_c3_transport_witness_constructor.json"


class TestP2156StrictD3C3TransportWitnessConstructor(unittest.TestCase):
    def test_export(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2156_s1106_strict_d3_c3_transport_witness_constructor.py")],
            check=True,
        )
        payload = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(payload["schema_version"], "p2156_s1106_v1")
        self.assertTrue(payload["gatekeeper_checks"]["witness_exported"])
        self.assertTrue(payload["strict_d3_c3_transport_witness_constructor"]["symbolic_zero_identity"])


if __name__ == "__main__":
    unittest.main()
