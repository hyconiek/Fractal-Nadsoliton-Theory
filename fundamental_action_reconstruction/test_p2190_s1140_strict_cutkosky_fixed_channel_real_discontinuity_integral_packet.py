from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2190S1140StrictCutkoskyFixedChannelRealDiscontinuityIntegralPacket(unittest.TestCase):
    def test_packet_and_numerics(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2190_s1140_strict_cutkosky_fixed_channel_real_discontinuity_integral_packet.py")],
            check=True,
        )
        data = json.loads(
            (
                G / "p2190_s1140_strict_cutkosky_fixed_channel_real_discontinuity_integral_packet.json"
            ).read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2190_s1140_v1")
        self.assertTrue(data["gatekeeper_checks"]["real_discontinuity_integral_exported"])
        self.assertIn("strict_cutkosky_fixed_channel_real_discontinuity_integral_packet", data)


if __name__ == "__main__":
    unittest.main()
