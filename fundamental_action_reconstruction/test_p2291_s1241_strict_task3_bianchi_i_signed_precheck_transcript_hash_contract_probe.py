from __future__ import annotations
import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2291S1241StrictTask3BianchiISignedPrecheckTranscriptHashContractProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2291_s1241_strict_task3_bianchi_i_signed_precheck_transcript_hash_contract_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2291_s1241_strict_task3_bianchi_i_signed_precheck_transcript_hash_contract_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2291_s1241_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["transcript_hash_exported"])
        self.assertTrue(g["hash_contract_exported"])
        self.assertTrue(g["decision_valid"])
        self.assertTrue(g["contract_status_valid"])


if __name__ == "__main__":
    unittest.main()
