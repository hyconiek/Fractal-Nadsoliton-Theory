from __future__ import annotations
import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2284S1234StrictTask3BianchiIQuantifiedPremiseTableAndTheoremBridgeDraftProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2284_s1234_strict_task3_bianchi_i_quantified_premise_table_and_theorem_bridge_draft_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2284_s1234_strict_task3_bianchi_i_quantified_premise_table_and_theorem_bridge_draft_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2284_s1234_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["premise_table_exported"])
        self.assertTrue(g["all_premise_ids_present"])
        self.assertTrue(g["draft_exported"])


if __name__ == "__main__":
    unittest.main()
