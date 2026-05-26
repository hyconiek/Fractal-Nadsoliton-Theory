from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2171P2172QW2191O1O3O4Progression(unittest.TestCase):
    def test_progression(self) -> None:
        subprocess.run([sys.executable, str(ROOT / "p2171_s1121_strict_qw2191_o1_o3_o4_dedicated_witness_packet.py")], check=True)
        subprocess.run([sys.executable, str(ROOT / "p2172_s1122_strict_qw2191_obligation_validator_o1_o3_o4_update.py")], check=True)

        d1 = json.loads((G / "p2171_s1121_strict_qw2191_o1_o3_o4_dedicated_witness_packet.json").read_text(encoding="utf-8"))
        d2 = json.loads((G / "p2172_s1122_strict_qw2191_obligation_validator_o1_o3_o4_update.json").read_text(encoding="utf-8"))

        self.assertEqual(d1["schema_version"], "p2171_s1121_v1")
        self.assertEqual(d2["schema_version"], "p2172_s1122_v1")
        self.assertTrue(d1["gatekeeper_checks"]["o1_pass"])
        self.assertTrue(d1["gatekeeper_checks"]["o3_pass"])
        self.assertTrue(d1["gatekeeper_checks"]["o4_pass"])
        self.assertTrue(d2["gatekeeper_checks"]["o1_promoted"])
        self.assertTrue(d2["gatekeeper_checks"]["o3_promoted"])
        self.assertTrue(d2["gatekeeper_checks"]["o4_promoted"])


if __name__ == "__main__":
    unittest.main()
