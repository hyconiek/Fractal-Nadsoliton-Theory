from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2204S1154StrictFrwBianchiSymbolicTransportOperatorMismatchCertificate(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2204_s1154_strict_frw_bianchi_symbolic_transport_operator_mismatch_certificate.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2204_s1154_strict_frw_bianchi_symbolic_transport_operator_mismatch_certificate.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2204_s1154_v1")
        self.assertTrue(data["gatekeeper_checks"]["symbolic_mismatch_certificate_exported"])
        self.assertTrue(data["gatekeeper_checks"]["symbolic_zero_mismatch_at_m_eq_1"])


if __name__ == "__main__":
    unittest.main()
