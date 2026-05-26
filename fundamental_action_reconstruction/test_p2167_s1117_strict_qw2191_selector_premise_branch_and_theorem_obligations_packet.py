from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p2167_s1117_strict_qw2191_selector_premise_branch_and_theorem_obligations_packet.json"


class TestP2167StrictQW2191SelectorPremiseBranchAndTheoremObligationsPacket(unittest.TestCase):
    def test_export(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2167_s1117_strict_qw2191_selector_premise_branch_and_theorem_obligations_packet.py")],
            check=True,
        )
        payload = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(payload["schema_version"], "p2167_s1117_v1")
        self.assertTrue(payload["gatekeeper_checks"]["branch_packet_exported"])
        self.assertTrue(payload["gatekeeper_checks"]["branch_instantiated"])
        self.assertGreater(len(payload["strict_qw2191_selector_premise_branch_and_theorem_obligations_packet"]["open_obligations"]), 0)


if __name__ == "__main__":
    unittest.main()
