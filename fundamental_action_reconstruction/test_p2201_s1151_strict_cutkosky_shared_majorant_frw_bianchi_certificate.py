from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2201S1151StrictCutkoskySharedMajorantFrwBianchiCertificate(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2201_s1151_strict_cutkosky_shared_majorant_frw_bianchi_certificate.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2201_s1151_strict_cutkosky_shared_majorant_frw_bianchi_certificate.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2201_s1151_v1")
        self.assertTrue(data["gatekeeper_checks"]["shared_majorant_certificate_exported"])


if __name__ == "__main__":
    unittest.main()
