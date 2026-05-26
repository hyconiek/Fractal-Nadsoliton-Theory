from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p2137_s1087_strict_cmp2_real_extension_schema_delivery_packet.json"


class TestP2137StrictCmp2RealExtensionSchemaDeliveryPacket(unittest.TestCase):
    def test_p2137_exports_schema_delivery_packet(self) -> None:
        subprocess.run([sys.executable, str(ROOT / "p2137_s1087_strict_cmp2_real_extension_schema_delivery_packet.py")], check=True)

        d = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(d["schema_version"], "p2137_s1087_v1")
        self.assertEqual(d["status"], "OPEN_PARTIAL_PROGRESS_WITH_TRACE")
        self.assertIn(d["result_kind"], {
            "PASS_STRICT_CMP2_REAL_EXTENSION_SCHEMA_DELIVERY_PACKET_WITH_TRACE",
            "OPEN_STRICT_CMP2_REAL_EXTENSION_SCHEMA_DELIVERY_PACKET_BLOCKED",
        })
        self.assertTrue(d["gatekeeper_checks"]["template_exported"])


if __name__ == "__main__":
    unittest.main()
