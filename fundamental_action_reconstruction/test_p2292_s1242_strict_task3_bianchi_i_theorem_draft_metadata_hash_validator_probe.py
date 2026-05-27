from __future__ import annotations
import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2292S1242StrictTask3BianchiITheoremDraftMetadataHashValidatorProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2292_s1242_strict_task3_bianchi_i_theorem_draft_metadata_hash_validator_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2292_s1242_strict_task3_bianchi_i_theorem_draft_metadata_hash_validator_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2292_s1242_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["validator_exported"])
        self.assertTrue(g["required_hash_length_ok"])
        self.assertTrue(g["checks_exported"])
        self.assertTrue(g["decision_valid"])


if __name__ == "__main__":
    unittest.main()
