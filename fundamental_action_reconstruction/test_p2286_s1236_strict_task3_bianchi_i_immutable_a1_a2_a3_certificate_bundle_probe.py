from __future__ import annotations
import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2286S1236StrictTask3BianchiIImmutableA1A2A3CertificateBundleProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2286_s1236_strict_task3_bianchi_i_immutable_a1_a2_a3_certificate_bundle_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2286_s1236_strict_task3_bianchi_i_immutable_a1_a2_a3_certificate_bundle_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2286_s1236_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["bundle_exported"])
        self.assertTrue(g["sha256_length_ok"])
        self.assertTrue(g["verifier_pass_boolean"])
        self.assertTrue(g["consistency_has_A1_A2_A3"])


if __name__ == "__main__":
    unittest.main()
