from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2191S1141StrictCutkoskyUncertaintyBandsAndUrResidualEnvelope(unittest.TestCase):
    def test_packet_and_envelope(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2191_s1141_strict_cutkosky_uncertainty_bands_and_ur_residual_envelope.py")],
            check=True,
        )
        data = json.loads(
            (
                G / "p2191_s1141_strict_cutkosky_uncertainty_bands_and_ur_residual_envelope.json"
            ).read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2191_s1141_v1")
        self.assertTrue(data["gatekeeper_checks"]["uncertainty_band_table_exported"])
        witness = data["strict_cutkosky_uncertainty_bands_and_ur_residual_envelope"]
        self.assertGreaterEqual(len(witness["uncertainty_band_table"]), 1)


if __name__ == "__main__":
    unittest.main()
