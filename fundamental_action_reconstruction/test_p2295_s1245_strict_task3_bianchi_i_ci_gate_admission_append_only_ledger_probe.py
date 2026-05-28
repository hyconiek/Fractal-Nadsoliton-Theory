from __future__ import annotations
import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2295S1245StrictTask3BianchiICiGateAdmissionAppendOnlyLedgerProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2295_s1245_strict_task3_bianchi_i_ci_gate_admission_append_only_ledger_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2295_s1245_strict_task3_bianchi_i_ci_gate_admission_append_only_ledger_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2295_s1245_v1")
        probe = data["strict_task3_bianchi_i_ci_gate_admission_append_only_ledger_probe"]
        self.assertTrue(probe["ledger_policy"]["append_only"])
        self.assertEqual(probe["ledger"][0]["run_index"], 1)
        self.assertEqual(probe["ledger"][0]["previous_entry_sha256"], None)
        self.assertEqual(probe["ledger"][0]["gate_decision"], "CI_GATE_BLOCK")
        self.assertEqual(probe["ledger"][0]["theorem_attempt_decision"], "THEOREM_DRAFT_HOLD")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["ledger_exported"])
        self.assertTrue(g["run_indices_monotonic"])
        self.assertTrue(g["first_run_index_is_one"])
        self.assertTrue(g["source_hashes_present"])
        self.assertTrue(g["entry_hash_length_ok"])
        self.assertTrue(g["ledger_fingerprint_length_ok"])
        self.assertTrue(g["mismatch_evidence_pointers_present"])
        self.assertTrue(g["gate_decision_preserved"])
        self.assertTrue(g["theorem_attempt_decision_preserved"])
        entry = dict(probe["ledger"][0])
        entry_hash = entry.pop("entry_sha256")
        recomputed = hashlib.sha256(
            json.dumps(entry, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
        ).hexdigest()
        self.assertEqual(entry_hash, recomputed)
        self.assertTrue(g["no_bridge_theorem_claimed"])
        self.assertTrue(g["no_selector_closure_claimed"])
        self.assertTrue(g["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
