from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2169P2170QW2191O2O5Progression(unittest.TestCase):
    def test_progression(self) -> None:
        subprocess.run([sys.executable, str(ROOT / "p2169_s1119_strict_qw2191_o2_o5_dedicated_witness_packet.py")], check=True)
        subprocess.run([sys.executable, str(ROOT / "p2170_s1120_strict_qw2191_obligation_validator_o2_o5_update.py")], check=True)

        d9 = json.loads((G / "p2169_s1119_strict_qw2191_o2_o5_dedicated_witness_packet.json").read_text(encoding="utf-8"))
        d0 = json.loads((G / "p2170_s1120_strict_qw2191_obligation_validator_o2_o5_update.json").read_text(encoding="utf-8"))

        self.assertEqual(d9["schema_version"], "p2169_s1119_v1")
        self.assertEqual(d0["schema_version"], "p2170_s1120_v1")
        self.assertTrue(d9["gatekeeper_checks"]["o2_pass"])
        self.assertTrue(d9["gatekeeper_checks"]["o5_pass"])
        self.assertTrue(d0["gatekeeper_checks"]["o2_promoted"])
        self.assertTrue(d0["gatekeeper_checks"]["o5_promoted"])


if __name__ == "__main__":
    unittest.main()
