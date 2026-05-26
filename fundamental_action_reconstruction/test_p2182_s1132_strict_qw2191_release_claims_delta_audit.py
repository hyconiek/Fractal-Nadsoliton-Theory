from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2182S1132StrictQW2191ReleaseClaimsDeltaAudit(unittest.TestCase):
    def test_release_delta_audit(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2182_s1132_strict_qw2191_release_claims_delta_audit.py")],
            check=True,
        )
        d = json.loads((G / "p2182_s1132_strict_qw2191_release_claims_delta_audit.json").read_text(encoding="utf-8"))
        self.assertEqual(d["schema_version"], "p2182_s1132_v1")
        self.assertTrue(d["gatekeeper_checks"]["release_delta_audit_exported"])
        self.assertIn("delta_findings", d["strict_qw2191_release_claims_delta_audit"])


if __name__ == "__main__":
    unittest.main()
