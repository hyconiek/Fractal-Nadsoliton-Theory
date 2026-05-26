from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2205S1155StrictNuBranchTransportVanishingConstraintsProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2205_s1155_strict_nu_branch_transport_vanishing_constraints_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2205_s1155_strict_nu_branch_transport_vanishing_constraints_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2205_s1155_v1")
        self.assertTrue(data["gatekeeper_checks"]["vanishing_constraints_probe_exported"])
        rows = data["strict_nu_branch_transport_vanishing_constraints_probe"]["numeric_local_scaling_rows"]
        self.assertGreaterEqual(len(rows), 5)


if __name__ == "__main__":
    unittest.main()
