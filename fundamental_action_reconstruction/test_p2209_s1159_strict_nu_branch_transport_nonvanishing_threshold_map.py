from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2209S1159StrictNuBranchTransportNonvanishingThresholdMap(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2209_s1159_strict_nu_branch_transport_nonvanishing_threshold_map.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2209_s1159_strict_nu_branch_transport_nonvanishing_threshold_map.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2209_s1159_v1")
        self.assertTrue(data["gatekeeper_checks"]["threshold_map_exported"])
        self.assertTrue(data["gatekeeper_checks"]["positive_linear_coeff"])
        self.assertGreater(len(data["strict_nu_branch_transport_nonvanishing_threshold_map"]["threshold_map_rows"]), 0)


if __name__ == "__main__":
    unittest.main()
