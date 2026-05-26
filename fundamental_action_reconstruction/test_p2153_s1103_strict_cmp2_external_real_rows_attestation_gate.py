from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p2153_s1103_strict_cmp2_external_real_rows_attestation_gate.json"


class TestP2153StrictCmp2ExternalRealRowsAttestationGate(unittest.TestCase):
    def test_export(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2153_s1103_strict_cmp2_external_real_rows_attestation_gate.py")],
            check=True,
        )
        payload = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(payload["schema_version"], "p2153_s1103_v1")
        self.assertTrue(payload["gatekeeper_checks"]["attestation_gate_exported"])


if __name__ == "__main__":
    unittest.main()
