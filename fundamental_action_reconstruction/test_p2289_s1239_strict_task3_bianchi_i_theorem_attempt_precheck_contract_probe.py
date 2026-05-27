from __future__ import annotations
import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2289S1239StrictTask3BianchiITheoremAttemptPrecheckContractProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2289_s1239_strict_task3_bianchi_i_theorem_attempt_precheck_contract_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2289_s1239_strict_task3_bianchi_i_theorem_attempt_precheck_contract_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2289_s1239_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["precheck_contract_exported"])
        self.assertTrue(g["fingerprint_length_ok"])
        self.assertTrue(g["observed_gating_decision_valid"])
        self.assertTrue(g["contract_status_valid"])


if __name__ == "__main__":
    unittest.main()
