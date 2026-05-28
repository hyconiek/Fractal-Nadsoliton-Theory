from __future__ import annotations
import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2288S1238StrictTask3BianchiICertificateChainIndexGatingLedgerProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2288_s1238_strict_task3_bianchi_i_certificate_chain_index_gating_ledger_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2288_s1238_strict_task3_bianchi_i_certificate_chain_index_gating_ledger_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2288_s1238_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["chain_index_exported"])
        self.assertTrue(g["bundle_hash_length_ok"])
        self.assertTrue(g["chain_fingerprint_length_ok"])
        self.assertTrue(g["gating_decision_present"])
        self.assertTrue(g["block_if_verifier_fails"])
        rec = data["strict_task3_bianchi_i_certificate_chain_index_gating_ledger_probe"]["chain_index_record"]
        self.assertFalse(rec["verifier_pass"])
        self.assertEqual(rec["gating_decision"], "BLOCK_THEOREM_ATTEMPT_PRECHECK")


if __name__ == "__main__":
    unittest.main()
