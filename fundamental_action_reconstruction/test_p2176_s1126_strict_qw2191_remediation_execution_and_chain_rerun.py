from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2176S1126StrictQW2191RemediationExecutionAndChainRerun(unittest.TestCase):
    def test_remediation_execution(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2176_s1126_strict_qw2191_remediation_execution_and_chain_rerun.py")],
            check=True,
        )
        d = json.loads(
            (G / "p2176_s1126_strict_qw2191_remediation_execution_and_chain_rerun.json").read_text(encoding="utf-8")
        )
        self.assertEqual(d["schema_version"], "p2176_s1126_v1")
        self.assertTrue(d["gatekeeper_checks"]["remediation_execution_exported"])
        summary = d["strict_qw2191_remediation_execution_and_chain_rerun"]["chain_rerun_summary"]
        self.assertIn("remaining_missing_count", summary)


if __name__ == "__main__":
    unittest.main()
