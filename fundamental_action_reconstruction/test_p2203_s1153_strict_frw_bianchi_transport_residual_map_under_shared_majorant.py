from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2203S1153StrictFrwBianchiTransportResidualMapUnderSharedMajorant(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2203_s1153_strict_frw_bianchi_transport_residual_map_under_shared_majorant.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2203_s1153_strict_frw_bianchi_transport_residual_map_under_shared_majorant.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2203_s1153_v1")
        self.assertTrue(data["gatekeeper_checks"]["transport_residual_map_exported"])
        rows = data["strict_frw_bianchi_transport_residual_map_under_shared_majorant"]["residual_map_rows"]
        self.assertGreaterEqual(len(rows), 10)


if __name__ == "__main__":
    unittest.main()
