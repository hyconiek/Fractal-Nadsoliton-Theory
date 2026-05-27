from __future__ import annotations
import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2287S1237StrictTask3BianchiICertificateBundleIntegrityRecheckProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2287_s1237_strict_task3_bianchi_i_certificate_bundle_integrity_recheck_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2287_s1237_strict_task3_bianchi_i_certificate_bundle_integrity_recheck_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2287_s1237_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["integrity_recheck_exported"])
        self.assertTrue(g["integrity_ok"])
        self.assertTrue(g["mutation_detected"])
        self.assertTrue(g["reported_hash_length_ok"])
        self.assertTrue(g["recomputed_hash_length_ok"])


if __name__ == "__main__":
    unittest.main()
