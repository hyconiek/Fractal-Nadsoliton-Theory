from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p2164_s1114_strict_qw2191_blocker_lane_baselined_entry_packet.json"


class TestP2164StrictQW2191BlockerLaneBaselinedEntryPacket(unittest.TestCase):
    def test_export(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2164_s1114_strict_qw2191_blocker_lane_baselined_entry_packet.py")],
            check=True,
        )
        payload = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(payload["schema_version"], "p2164_s1114_v1")
        self.assertTrue(payload["gatekeeper_checks"]["baselined_entry_exported"])
        self.assertEqual(
            payload["strict_qw2191_blocker_lane_baselined_entry_packet"]["selected_next_blocker_lane"],
            "QW-2191_selector_obstruction",
        )


if __name__ == "__main__":
    unittest.main()
