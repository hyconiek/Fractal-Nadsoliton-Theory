from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2198S1148StrictReleaseGatePolicyContractAndContradictionDeltaBridge(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [
                sys.executable,
                str(ROOT / "p2198_s1148_strict_release_gate_policy_contract_and_contradiction_delta_bridge.py"),
            ],
            check=True,
        )
        data = json.loads(
            (
                G
                / "p2198_s1148_strict_release_gate_policy_contract_and_contradiction_delta_bridge.json"
            ).read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2198_s1148_v1")
        self.assertTrue(data["gatekeeper_checks"]["policy_contract_delta_bridge_exported"])


if __name__ == "__main__":
    unittest.main()
