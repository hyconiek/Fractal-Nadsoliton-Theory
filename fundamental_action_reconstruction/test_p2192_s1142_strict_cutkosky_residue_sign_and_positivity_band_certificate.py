from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2192S1142StrictCutkoskyResidueSignAndPositivityBandCertificate(unittest.TestCase):
    def test_packet_and_certificate(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2192_s1142_strict_cutkosky_residue_sign_and_positivity_band_certificate.py")],
            check=True,
        )
        data = json.loads(
            (
                G / "p2192_s1142_strict_cutkosky_residue_sign_and_positivity_band_certificate.json"
            ).read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2192_s1142_v1")
        self.assertTrue(data["gatekeeper_checks"]["certificate_exported"])
        cert = data["strict_cutkosky_residue_sign_and_positivity_band_certificate"]
        self.assertGreater(cert["sample_count"], 0)


if __name__ == "__main__":
    unittest.main()
