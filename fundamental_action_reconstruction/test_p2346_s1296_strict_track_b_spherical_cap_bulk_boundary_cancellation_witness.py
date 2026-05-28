from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2346_s1296_strict_track_b_spherical_cap_bulk_boundary_cancellation_witness.py"
OUT = ROOT / "generated" / "p2346_s1296_strict_track_b_spherical_cap_bulk_boundary_cancellation_witness.json"
MD = ROOT / "generated" / "p2346_s1296_strict_track_b_spherical_cap_bulk_boundary_cancellation_witness.md"


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2346SphericalCapBulkBoundaryCancellationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["strict_track_b_spherical_cap_bulk_boundary_cancellation_witness_probe"]
        cls.witness = cls.probe["track_B_spherical_cap_bulk_boundary_cancellation_witness"]
        cls.cap = cls.witness["spherical_cap_class"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2346_s1296_v1")
        self.assertEqual(self.payload["packet_id"], "P2346")
        self.assertEqual(self.payload["stage_id"], "S1296")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_spherical_cap_bulk_boundary_values(self) -> None:
        self.assertEqual(self.witness["ledger_coefficient_b_GB"], "13152087*log(2)/(320000000*pi**2)")
        self.assertEqual(self.witness["per_euler_pairing"], "13152087*log(2)/10000000")
        self.assertEqual(self.cap["unit_s4_gb_density"], "24")
        self.assertEqual(self.cap["cap_volume"], "2*pi**2*(c - 1)**2*(c + 2)/3")
        self.assertEqual(self.cap["bulk_gb_integral"], "16*pi**2*(c - 1)**2*(c + 2)")
        self.assertEqual(self.cap["boundary_transgression"], "-16*pi**2*c*(c**2 - 3)")
        self.assertEqual(self.cap["bulk_plus_boundary"], "32*pi**2")
        self.assertEqual(self.cap["target_topological_number_32pi2_chi"], "32*pi**2")
        self.assertEqual(self.cap["cancellation_residual"], "0")
        self.assertEqual(self.cap["c_derivative_residual"], "0")
        self.assertEqual(self.cap["flat_limit_bulk_c1"], "0")
        self.assertEqual(self.cap["flat_limit_boundary_c1"], "32*pi**2")
        self.assertEqual(self.cap["hemisphere_bulk_c0"], "32*pi**2")
        self.assertEqual(self.cap["hemisphere_boundary_c0"], "0")
        self.assertEqual(self.cap["sample_c_minus_half_total"], "32*pi**2")
        self.assertEqual(self.cap["normalized_track_B_pairing"], "13152087*log(2)/10000000")
        self.assertEqual(self.cap["normalized_pairing_residual"], "0")

    def test_dependencies_gatekeepers_and_fingerprint(self) -> None:
        deps = self.probe["current_export_dependencies"]
        self.assertTrue(deps["p2341_flat_limit_boundary_rederived"])
        self.assertTrue(deps["p2345_pairing_rederived"])
        for key, value in self.payload["gatekeeper_checks"].items():
            self.assertTrue(value, key)
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("general curved bulk-plus-boundary theorem for arbitrary four-manifolds", theorem["not_licensed"])
        self.assertIn("selector premise or QW-2191 selector discharge", theorem["not_licensed"])
        self.assertIn("ToE closure", theorem["not_licensed"])


if __name__ == "__main__":
    unittest.main()
