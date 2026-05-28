from __future__ import annotations
import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2285S1235StrictTask3BianchiIMachineCheckableA1A2A3VerifierProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2285_s1235_strict_task3_bianchi_i_machine_checkable_a1_a2_a3_verifier_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2285_s1235_strict_task3_bianchi_i_machine_checkable_a1_a2_a3_verifier_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2285_s1235_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["consistency_rows_exported"])
        self.assertTrue(g["all_consistency_flags_boolean"])
        self.assertTrue(g["contains_verifier_pass"])
        probe = data["strict_task3_bianchi_i_machine_checkable_a1_a2_a3_verifier_probe"]
        self.assertEqual(probe["recomputed"], {"A1": False, "A2": False, "A3": False})
        self.assertFalse(probe["verifier_pass"])


if __name__ == "__main__":
    unittest.main()
