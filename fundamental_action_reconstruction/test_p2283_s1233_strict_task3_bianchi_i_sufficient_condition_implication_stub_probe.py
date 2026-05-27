from __future__ import annotations
import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2283S1233StrictTask3BianchiISufficientConditionImplicationStubProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2283_s1233_strict_task3_bianchi_i_sufficient_condition_implication_stub_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2283_s1233_strict_task3_bianchi_i_sufficient_condition_implication_stub_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2283_s1233_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["assumption_rows_exported"])
        self.assertTrue(g["contains_A1"])
        self.assertTrue(g["contains_A2"])
        self.assertTrue(g["contains_A3"])


if __name__ == "__main__":
    unittest.main()
