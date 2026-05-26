from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2188S1138StrictQW2191ReleaseGateContractFromCiPolicy(unittest.TestCase):
    def test_release_gate_contract_export(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2188_s1138_strict_qw2191_release_gate_contract_from_ci_policy.py")],
            check=True,
        )
        data = json.loads(
            (
                G / "p2188_s1138_strict_qw2191_release_gate_contract_from_ci_policy.json"
            ).read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2188_s1138_v1")
        self.assertTrue(data["gatekeeper_checks"]["release_gate_contract_exported"])
        contract = data["strict_qw2191_release_gate_contract_from_ci_policy"]
        self.assertIn(contract["output"]["release_gate_decision"], {"RELEASE_GATE_OPEN", "RELEASE_GATE_OPEN_WITH_NOTICE", "RELEASE_GATE_BLOCKED"})


if __name__ == "__main__":
    unittest.main()
